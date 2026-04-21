estimate_unpenalized_FCIR <- function(Y, X, Tr){
  # Y: N x P binary response matrix
  # X: N x L environment matrix (first column should be 1s for intercept)
  # Tr: P x K species traits matrix
  
  N = nrow(Y)
  P = ncol(Y)
  L = ncol(X)
  K = ncol(Tr)
  
  # 1. Pre-compute Trait Differences
  # Calculate the absolute trait differences Delta_{jj'} for all species pairs
  Delta = array(0, dim = c(P, P, K))
  for(j in 1:P) {
    for(j_prime in 1:P) {
      Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
    }
  }
  
  # 2. Initialize pseudo-likelihood response vector and design matrix
  n_obs = N * P               # Total number of observations (rows) for pseudo-likelihood estimation
  n_params = 2*L + 2*L*K      # Total number of parameters in the full model: beta_0, vec(B), alpha_0, vec(A)
  
  glm_Y = numeric(n_obs)
  glm_X = matrix(0, nrow = n_obs, ncol = n_params)
  
  # 3. Construct the designmatrix X (corresponding to the Kronecker products and block structure in the paper)
  row_idx = 1
  for(s in 1:N){
    x_s = X[s, ]     # Environment vector for current site s
    y_s = Y[s, ]     # Distribution of all species at current site s
    R_s = sum(y_s)   # Species richness at current site s (used for neighbor summation)
    
    for(j in 1:P){
      # Response variable: presence of species j at site s
      glm_Y[row_idx] = y_s[j]
      
      # --- Main Effects ---
      t_j = Tr[j, ]
      comp1_beta0 = x_s
      comp2_B     = kronecker(t_j, x_s) # Equivalent to t_j \otimes x_s
      
      # --- Interaction Effects ---
      neighbor_sum = R_s - y_s[j]       # Equivalent to sum_{j' \neq j} y_sj'
      comp3_alpha0 = neighbor_sum * x_s
      
      # Calculate the neighbor trait difference weighted sum (w_sj) for focal species j
      w_sj = numeric(K)
      for(j_prime in 1:P){
        if(j_prime != j && y_s[j_prime] == 1){
          w_sj = w_sj + Delta[j, j_prime, ]
        }
      }
      comp4_A = kronecker(w_sj, x_s)    # Equivalent to w_sj \otimes x_s
      
      # Concatenate the four blocks into one row for the current observation
      glm_X[row_idx, ] = c(comp1_beta0, comp2_B, comp3_alpha0, comp4_A)
      
      row_idx = row_idx + 1
    }
  }
  
  # 4. Fit Unpenalized Logistic Regression
  # Note: '+ 0' is added to the formula because the intercept is already implicitly 
  # included in the design matrix (first column of x_s is 1)
  logistic_reg = glm(glm_Y ~ glm_X + 0, family = binomial, control = glm.control(maxit = 100))
  
  # 5. Extract and reshape parameters
  est_coefs = logistic_reg$coefficients
  
  idx = 1
  hat_beta_0  = est_coefs[idx:(idx + L - 1)]; idx = idx + L
  hat_B_vec   = est_coefs[idx:(idx + L*K - 1)]; idx = idx + L*K
  hat_alpha_0 = est_coefs[idx:(idx + L - 1)]; idx = idx + L
  hat_A_vec   = est_coefs[idx:(idx + L*K - 1)]
  
  # Reshape vectors back into L x K fourth-corner matrices
  hat_B_mat = matrix(hat_B_vec, nrow = L, ncol = K)
  hat_A_mat = matrix(hat_A_vec, nrow = L, ncol = K)
  
  return(list(
    beta_0 = hat_beta_0,
    B_mat = hat_B_mat,
    alpha_0 = hat_alpha_0,
    A_mat = hat_A_mat,
    glm_model = logistic_reg # Return the original model to allow checking p-values or summary
  ))
}