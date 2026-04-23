library(IsingSampler)

generate_dense_fcir_data <- function(N, P, L, K, B_reps, seed, filename){
  # This function generates site-varying Ising data based on the FCIR formulation
  
  set.seed(seed)
  
  # 1. Initialize Matrices
  Y = array(data = NA, dim = c(N, P, B_reps))
  X = matrix(rnorm(N * L), nrow = N, ncol = L)
  X[,1] = 1 
  Tr = matrix(runif(P * K, min = -1, max = 1), nrow = P, ncol = K)
  
  # 2. Helper function to generate sparse parameters
  generate_sparse_params <- function(n_elements, prob_zero, min_mag, max_mag) {
    is_zero = rbinom(n_elements, 1, prob_zero)
    mags = runif(n_elements, min_mag, max_mag)
    signs = sample(c(-1, 1), n_elements, replace = TRUE)
    return((1 - is_zero) * signs * mags)
  }
  
  # 3. Generate Main Effect Parameters
  beta_0 = generate_sparse_params(L, prob_zero = 0.3, min_mag = 0.5, max_mag = 1.5)
  B_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.3, min_mag = 0.2, max_mag = 0.5), nrow = L, ncol = K)
  
  # 4. Generate Interaction Effect Parameters (NO Adj MATRIX)
  alpha_0 = generate_sparse_params(L, prob_zero = 0.3, min_mag = 0.4, max_mag = 1.0)
  A_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.5, min_mag = 0.3, max_mag = 0.8), nrow = L, ncol = K)
  
  # 5. Pre-compute pairwise trait differences
  Delta = array(0, dim = c(P, P, K))
  for(j in 1:P) {
    for(j_prime in 1:P) {
      Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
    }
  }
  
  # 6. MCMC Sampling for each replication
  Beta_temp = 1 
  
  for(b in 1:B_reps){
    Y_b = matrix(NA, nrow = N, ncol = P)
    
    for(s in 1:N){
      x_s = X[s, ] 
      
      # Step A: Calculate site-specific main effects (theta_jj)
      theta_jj_s = numeric(P)
      for(j in 1:P){
        t_j = Tr[j, ]
        theta_jj_s[j] = sum(x_s * beta_0) + sum(x_s * (B_mat %*% t_j))
      }
      
      # Step B: Calculate site-specific interaction network (Theta)
      Theta_s = matrix(0, nrow = P, ncol = P)
      for(j in 1:(P-1)){
        for(j_prime in (j+1):P){
          delta_jj = Delta[j, j_prime, ]
          # Calculate edge value directly for ALL pairs
          edge_val = sum(x_s * alpha_0) + sum(x_s * (A_mat %*% delta_jj))
          Theta_s[j, j_prime] = edge_val
          Theta_s[j_prime, j] = edge_val
        }
      }
      
      # Step C: Sample response for this specific site
      sampled_y = IsingSampler(1, Theta_s, theta_jj_s, Beta_temp, 1000/P, 
                               responses = c(0L, 1L), method = "MH")
      Y_b[s, ] = sampled_y
    }
    
    Y[,,b] = Y_b
  }
  
  # 7. Save output 
  save(Y, X, Tr, beta_0, B_mat, alpha_0, A_mat, 
       N, P, L, K, B_reps, seed, file = filename)
  
}
