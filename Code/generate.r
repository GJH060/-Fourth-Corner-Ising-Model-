library(IsingSampler)

generate_fcir_data <- function(N, P, L, K, B_reps, seed, p_edge, filename){
  # This function generates site-varying Ising data based on the Fourth-Corner Ising Regression (FCIR) formulation.
  # It integrates environment covariates (X) and species functional traits (T) to dynamically shape
  # both the main occurrence effects and the pairwise species interaction network.
  
  set.seed(seed)
  
  # 1. Initialize Matrices
  Y = array(data = NA, dim = c(N, P, B_reps))
  
  # Generate Environment Covariates (X): N x L
  # Includes intercept (col 1) + independent standard normal variables
  X = matrix(rnorm(N * L), nrow = N, ncol = L)
  X[,1] = 1 
  
  # Generate Species Traits (Tr): P x K
  # Drawn independently from Uniform(-1, 1)
  Tr = matrix(runif(P * K, min = -1, max = 1), nrow = P, ncol = K)
  
  # 2. Helper function to generate sparse parameters with specific magnitude bounds
  generate_sparse_params <- function(n_elements, prob_zero, min_mag, max_mag) {
    is_zero = rbinom(n_elements, 1, prob_zero)
    mags = runif(n_elements, min_mag, max_mag)
    signs = sample(c(-1, 1), n_elements, replace = TRUE)
    return((1 - is_zero) * signs * mags)
  }
  
  # 3. Generate Main Effect Parameters
  # beta_0: L-dimensional vector for environment main effects, U(0.5, 1.5)
  beta_0 = generate_sparse_params(L, prob_zero = 0.3, min_mag = 0.5, max_mag = 1.5)
  
  # B_mat: L x K matrix for trait-environment main effects (Fourth-corner B), U(0.2, 0.5)
  B_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.3, min_mag = 0.2, max_mag = 0.5), nrow = L, ncol = K)
  
  # alpha_0: L-dimensional vector for base environment interaction, U(0.4, 1.0)
  alpha_0 = generate_sparse_params(L, prob_zero = 0.3, min_mag = 0.4, max_mag = 1.0)
  
  # A_mat: L x K matrix for trait-environment interaction (Fourth-corner A), U(0.3, 0.8)
  # Applying stricter sparsity here as defined in the simulation text
  A_mat = matrix(generate_sparse_params(L * K, prob_zero = 0.5, min_mag = 0.3, max_mag = 0.8), nrow = L, ncol = K)
  
  # 4. Pre-compute pairwise trait differences (element-wise absolute difference)
  # Delta array of dimension P x P x K
  Delta = array(0, dim = c(P, P, K))
  for(j in 1:P) {
    for(j_prime in 1:P) {
      Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
    }
  }
  
  # 5. Generate response data for each replication
  Beta_temp = 1 # Inverse temperature
  
  for(b in 1:B_reps){
    Y_b = matrix(NA, nrow = N, ncol = P)
    
    # Due to site-varying covariates, we MUST calculate potentials and sample for each site individually
    for(s in 1:N){
      x_s = X[s, ] # Current site's environment vector (L x 1)
      
      # Step A: Calculate site-specific main effects (theta_jj)
      # Equation: theta_{s,jj} = x_s^T * beta_0 + x_s^T * B * t_j
      theta_jj_s = numeric(P)
      for(j in 1:P){
        t_j = Tr[j, ]
        theta_jj_s[j] = sum(x_s * beta_0) + sum(x_s * (B_mat %*% t_j))
      }
      
      # Step B: Calculate site-specific interaction network (Theta)
      # Equation: theta_{s,jj'} = x_s^T * alpha_0 + x_s^T * A * Delta_{jj'}
      Theta_s = matrix(0, nrow = P, ncol = P)
      for(j in 1:(P-1)){
        for(j_prime in (j+1):P){
          delta_jj = Delta[j, j_prime, ]
          edge_val = sum(x_s * alpha_0) + sum(x_s * (A_mat %*% delta_jj))
          Theta_s[j, j_prime] = edge_val
          Theta_s[j_prime, j] = edge_val
        }
      }
    }
      
      # Step C: Sample response for this specific site
      # n=1 because the specific Theta_s and theta_jj_s are unique to site s
      sampled_y = IsingSampler(1, Theta_s, theta_jj_s, Beta_temp, 1000/P, 
                               responses = c(0L, 1L), method = "MH")
      Y_b[s, ] = sampled_y
    }
    
    Y[,,b] = Y_b
  }