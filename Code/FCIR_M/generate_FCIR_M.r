library(IsingSampler)

generate_fcir_M_data <- function(N, P, L, K, B_reps, seed, filename){
  # Restricted Model: FCIR_M
  # Traits and environment drive the main effects only.
  # Interaction network is simplified to static, site-independent residual associations.
  
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
  
  # Main Effect Parameters (Retained from full model)
  beta_0 = generate_sparse_params(L, prob_zero = 0.3, min_mag = 0.5, max_mag = 1.5)
  B_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.3, min_mag = 0.2, max_mag = 0.5), nrow = L, ncol = K)
  
  # --- MODIFICATION 1: Static Residual Interaction Network ---
  # Theta_int: P x P symmetric matrix independent of environment and traits.
  Theta_int = matrix(0, nrow = P, ncol = P)
  n_edges = P * (P - 1) / 2
  # Generate sparse static edge potentials 
  upper_tri_vals = generate_sparse_params(n_edges, prob_zero = 0.6, min_mag = 0.4, max_mag = 1.0)
  Theta_int[upper.tri(Theta_int)] = upper_tri_vals
  Theta_int[lower.tri(Theta_int)] = t(Theta_int)[lower.tri(Theta_int)]
  
  # (No Delta computation needed for FCIR_M)
  
  Beta_temp = 1 
  
  for(b in 1:B_reps){
    Y_b = matrix(NA, nrow = N, ncol = P)
    for(s in 1:N){
      x_s = X[s, ] 
      
      # Step A: Calculate site-specific main effects (Dynamic)
      theta_jj_s = numeric(P)
      for(j in 1:P){
        t_j = Tr[j, ]
        theta_jj_s[j] = sum(x_s * beta_0) + sum(x_s * (B_mat %*% t_j))
      }
      
      # --- MODIFICATION 2: Apply static network ---
      # Equation: theta_{s,jj'} = theta_{jj'}
      Theta_s = Theta_int 
      
      # Step C: Sample response 
      sampled_y = IsingSampler(1, Theta_s, theta_jj_s, Beta_temp, 1000/P, 
                               responses = c(0L, 1L), method = "MH")
      Y_b[s, ] = sampled_y
    }
    Y[,,b] = Y_b
  }
  save(Y, X, Tr, beta_0, B_mat, Theta_int, N, P, L, K, B_reps, seed, file = filename)
}