estimate_unpenalized_FCIR_I <- function(Y, X, Tr){
  # Y: N x P binary response matrix
  # X: N x L environment matrix (first column is 1s for intercept)
  # Tr: P x K species traits matrix
  
  N = nrow(Y)
  P = ncol(Y)
  L = ncol(X)
  K = ncol(Tr)
  
  # 1. Pre-compute Trait Differences (Needed for Interaction)
  Delta = array(0, dim = c(P, P, K))
  for(j in 1:P) {
    for(j_prime in 1:P) {
      Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
    }
  }
  
  # 2. Initialize design matrix for FCIR_I
  n_obs = N * P
  # Parameters: P species-specific beta vectors (P*L), alpha_0 (L), vec(A) (L*K)
  n_params = P*L + L + L*K      
  
  glm_Y = numeric(n_obs)
  glm_X = matrix(0, nrow = n_obs, ncol = n_params)
  
  row_idx = 1
  for(s in 1:N){
    x_s = X[s, ]
    y_s = Y[s, ]
    R_s = sum(y_s)
    
    for(j in 1:P){
      glm_Y[row_idx] = y_s[j]
      
      # --- MODIFICATION 1: Species-specific Main Effects ---
      # Create a sparse vector where only the L columns for species j contain x_s
      comp1_beta_species = numeric(P * L)
      start_col = (j - 1) * L + 1
      end_col = j * L
      comp1_beta_species[start_col:end_col] = x_s
      
      # --- Interaction Effects (Same as full model) ---
      neighbor_sum = R_s - y_s[j]
      comp2_alpha0 = neighbor_sum * x_s
      
      w_sj = numeric(K)
      for(j_prime in 1:P){
        if(j_prime != j && y_s[j_prime] == 1){
          w_sj = w_sj + Delta[j, j_prime, ]
        }
      }
      comp3_A = kronecker(w_sj, x_s)
      
      # Concatenate blocks
      glm_X[row_idx, ] = c(comp1_beta_species, comp2_alpha0, comp3_A)
      
      row_idx = row_idx + 1
    }
  }
  
  # 3. Fit Unpenalized Logistic Regression
  logistic_reg = glm(glm_Y ~ glm_X + 0, family = binomial, control = glm.control(maxit = 100))
  est_coefs = logistic_reg$coefficients
  
  # 4. Extract and reshape parameters
  idx = 1
  hat_Beta_vec = est_coefs[idx:(idx + P*L - 1)]; idx = idx + P*L
  hat_alpha_0  = est_coefs[idx:(idx + L - 1)];   idx = idx + L
  hat_A_vec    = est_coefs[idx:(idx + L*K - 1)]
  
  # Reshape Beta vector into a P x L matrix (Row j is beta_j)
  # R fills matrices column-wise, so we create L x P first, then transpose
  hat_Beta_mat = t(matrix(hat_Beta_vec, nrow = L, ncol = P))
  hat_A_mat = matrix(hat_A_vec, nrow = L, ncol = K)
  
  return(list(
    Beta_mat = hat_Beta_mat,
    alpha_0 = hat_alpha_0,
    A_mat = hat_A_mat,
    glm_model = logistic_reg
  ))
}