library(IsingSampler)

generate_dense_fcir_data <- function(N, P, L, K, B_reps, seed, filename){
  # This function generates site-varying Ising data based on the FCIR formulation
  
  set.seed(seed)
  
  # 2. Helper function to generate sparse parameters
  generate_sparse_params <- function(n_elements, prob_zero, min_mag, max_mag) {
    is_zero = rbinom(n_elements, 1, prob_zero)
    mags = runif(n_elements, min_mag, max_mag)
    signs = sample(c(-1, 1), n_elements, replace = TRUE)
    return((1 - is_zero) * signs * mags)
  }
  
  # 3. Generate Main Effect Parameters
  beta_0 = generate_sparse_params(L, prob_zero = 0, min_mag = 0.5, max_mag = 1.5)
  B_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.3, min_mag = 0.3, max_mag = 0.8), nrow = L, ncol = K)
  
  # 4. Generate Interaction Effect Parameters (NO Adj MATRIX)
  alpha_0 = generate_sparse_params(L, prob_zero = 0, min_mag = 0.4, max_mag = 1.0)*0.25
  A_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.5, min_mag = 0.3, max_mag = 0.8), nrow = L, ncol = K)*0.25
  
  
  # 1. Initialize Matrices
  Y = array(data = NA, dim = c(N, P, B_reps))
  X = matrix(rnorm(N * L), nrow = N, ncol = L)
  X[,1] = 1 
  Tr = matrix(runif(P * K, min = -1, max = 1), nrow = P, ncol = K)
  
  
  # 5. Pre-compute pairwise trait differences (fixed across replicates)
  Delta = array(0, dim = c(P, P, K))
  for(j in 1:P) {
    for(jp in 1:P) {
      Delta[j, jp, ] = abs(Tr[j, ] - Tr[jp, ])
    }
  }
  
  Beta_temp = 1 
  
  # 6. Pre-compute the site-specific fields once per site.
  # theta_jj_s (main-effect field) and Theta_s (interaction network) depend only
  # on X, Tr and the fixed parameters, NOT on the replicate index b, so building
  # them inside the b loop repeats the same O(P^2) work B_reps times. We compute
  # each site once and vectorise the pairwise construction:
  #   theta_jj_s[j]   = x_s . beta_0 + x_s . (B_mat %*% Tr[j, ])
  #   Theta_s[j, jp]  = x_s . alpha_0 + sum_k (x_s^T A_mat)[k] * Delta[j, jp, k]
  theta_jj_all = matrix(0, nrow = N, ncol = P)   # site s -> length-P field
  Theta_all = vector("list", N)                  # site s -> P x P interaction matrix
  for (s in 1:N) {
    x_s = X[s, ]
    theta_jj_all[s, ] = as.numeric(sum(x_s * beta_0) + x_s %*% B_mat %*% t(Tr))
    
    v_s = as.numeric(x_s %*% A_mat)              # length K: (x_s^T A_mat)
    Theta_s = matrix(sum(x_s * alpha_0), nrow = P, ncol = P)
    for (k in 1:K) Theta_s = Theta_s + Delta[, , k] * v_s[k]
    diag(Theta_s) = 0                            # no self-interaction
    Theta_all[[s]] = Theta_s
  }
  
  # 7. MCMC Sampling for each replication (the only b-dependent step)
  for(b in 1:B_reps){
    Y_b = matrix(NA, nrow = N, ncol = P)
    
    for(s in 1:N){
      # Sample response for this specific site
      Y_b[s, ] = IsingSampler(1, Theta_all[[s]], theta_jj_all[s, ], Beta_temp, 1000/P, 
                              responses = c(0L, 1L), method = "MH")
    }
    
    Y[,,b] = Y_b
  }
  
  # 8. Save output 
  save(Y, X, Tr, beta_0, B_mat, alpha_0, A_mat, 
       N, P, L, K, B_reps, seed, file = filename)
  
}
