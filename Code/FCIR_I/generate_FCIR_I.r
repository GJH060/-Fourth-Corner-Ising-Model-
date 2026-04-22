library(IsingSampler)

generate_fcir_I_data <- function(N, P, L, K, B_reps, seed, filename){
  # Restricted Model: FCIR_I 
  # Traits drive the species interaction network only. 
  # Main effects are unconstrained species-specific environmental responses.
  
  set.seed(seed)
  Y = array(data = NA, dim = c(N, P, B_reps))
  X = matrix(rnorm(N * L), nrow = N, ncol = L)
  X[,1] = 1 
  Tr = matrix(runif(P * K, min = -1, max = 1), nrow = P, ncol = K)
  
  generate_sparse_params <- function(n_elements, prob_zero, min_mag, max_mag) {
    is_zero = rbinom(n_elements, 1, prob_zero)
    mags = runif(n_elements, min_mag, max_mag)
    signs = sample(c(-1, 1), n_elements, replace = TRUE)
    return((1 - is_zero) * signs * mags)
  }
  
  # --- MODIFICATION 1: Species-specific Main Effects ---
  # Beta_mat: P x L matrix. Each row is beta_j for species j.
  # Parameter space expands from L(K+1) to P*L.
  Beta_mat = matrix(generate_sparse_params(P * L, prob_zero = 0.3, min_mag = 0.5, max_mag = 1.5), nrow = P, ncol = L)
  
  # Interaction Effect Parameters (Retained from full model)
  alpha_0 = generate_sparse_params(L, prob_zero = 0.3, min_mag = 0.4, max_mag = 1.0)
  A_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.5, min_mag = 0.3, max_mag = 0.8), nrow = L, ncol = K)
  
  Delta = array(0, dim = c(P, P, K))
  for(j in 1:P) {
    for(j_prime in 1:P) {
      Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
    }
  }
  
  Beta_temp = 1 
  
  for(b in 1:B_reps){
    Y_b = matrix(NA, nrow = N, ncol = P)
    for(s in 1:N){
      x_s = X[s, ] 
      
      # --- MODIFICATION 2: Calculate unconstrained main effects ---
      # Equation: theta_{s,jj} = x_s^T * beta_j
      theta_jj_s = numeric(P)
      for(j in 1:P){
        theta_jj_s[j] = sum(x_s * Beta_mat[j, ])
      }
      
      # Step B: Calculate site-specific interaction network (Dynamic)
      Theta_s = matrix(0, nrow = P, ncol = P)
      for(j in 1:(P-1)){
        for(j_prime in (j+1):P){
          delta_jj = Delta[j, j_prime, ]
          edge_val = sum(x_s * alpha_0) + sum(x_s * (A_mat %*% delta_jj))
          Theta_s[j, j_prime] = edge_val
          Theta_s[j_prime, j] = edge_val
        }
      }
      
      # Step C: Sample response (Fixed brace position)
      sampled_y = IsingSampler(1, Theta_s, theta_jj_s, Beta_temp, 1000/P, 
                               responses = c(0L, 1L), method = "MH")
      Y_b[s, ] = sampled_y
    }
    Y[,,b] = Y_b
  }
  save(Y, X, Tr, Beta_mat, alpha_0, A_mat, N, P, L, K, B_reps, seed, file = filename)
}