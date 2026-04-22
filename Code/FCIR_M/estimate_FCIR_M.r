estimate_unpenalized_FCIR_M <- function(Y, X, Tr){
  # Y: N x P binary response matrix
  # X: N x L environment matrix (first column is 1s for intercept)
  # Tr: P x K species traits matrix
  
  N = nrow(Y)
  P = ncol(Y)
  L = ncol(X)
  K = ncol(Tr)
  n_edges = P * (P - 1) / 2
  
  # 1. Pre-compute an Edge Index Matrix 
  # This mapping matrix quickly retrieves the 1D index (1 to n_edges) for any pair (j, j')
  pair_idx_mat = matrix(0, nrow = P, ncol = P)
  pair_idx_mat[upper.tri(pair_idx_mat)] = 1:n_edges
  pair_idx_mat[lower.tri(pair_idx_mat)] = t(pair_idx_mat)[lower.tri(pair_idx_mat)]
  
  # (No Delta computation needed for FCIR_M)
  
  # 2. Initialize design matrix for FCIR_M
  n_obs = N * P
  # Parameters: beta_0 (L), vec(B) (L*K), theta_int (n_edges)
  n_params = L + L*K + n_edges      
  
  glm_Y = numeric(n_obs)
  glm_X = matrix(0, nrow = n_obs, ncol = n_params)
  
  row_idx = 1
  for(s in 1:N){
    x_s = X[s, ]
    y_s = Y[s, ]
    
    for(j in 1:P){
      glm_Y[row_idx] = y_s[j]
      
      # --- Main Effects (Same as full model) ---
      t_j = Tr[j, ]
      comp1_beta0 = x_s
      comp2_B     = kronecker(t_j, x_s)
      
      # --- MODIFICATION 2: Static Edge-Incidence Vector ---
      # An indicator vector of length n_edges. Element is 1 if neighbor j' is present.
      comp3_theta_int = numeric(n_edges)
      for(j_prime in 1:P){
        if(j_prime != j && y_s[j_prime] == 1){
          edge_idx = pair_idx_mat[j, j_prime]
          comp3_theta_int[edge_idx] = 1 # Activate the corresponding edge parameter
        }
      }
      
      # Concatenate blocks
      glm_X[row_idx, ] = c(comp1_beta0, comp2_B, comp3_theta_int)
      
      row_idx = row_idx + 1
    }
  }
  
  # 3. Fit Unpenalized Logistic Regression
  logistic_reg = glm(glm_Y ~ glm_X + 0, family = binomial, control = glm.control(maxit = 100))
  est_coefs = logistic_reg$coefficients
  
  # 4. Extract and reshape parameters
  idx = 1
  hat_beta_0 = est_coefs[idx:(idx + L - 1)]; idx = idx + L
  hat_B_vec  = est_coefs[idx:(idx + L*K - 1)]; idx = idx + L*K
  hat_theta_int_vec = est_coefs[idx:(idx + n_edges - 1)]
  
  hat_B_mat = matrix(hat_B_vec, nrow = L, ncol = K)
  
  # Reconstruct the Theta_int matrix from the estimated 1D vector
  hat_Theta_int = matrix(0, nrow = P, ncol = P)
  hat_Theta_int[upper.tri(hat_Theta_int)] = hat_theta_int_vec
  hat_Theta_int[lower.tri(hat_Theta_int)] = t(hat_Theta_int)[lower.tri(hat_Theta_int)]
  
  return(list(
    beta_0 = hat_beta_0,
    B_mat = hat_B_mat,
    Theta_int = hat_Theta_int,
    glm_model = logistic_reg
  ))
}