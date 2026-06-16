library(glmnet)

estimate_penalized_FCIR <- function(Y, X, Tr, alpha = 1, lambda = "lambda.min", use_cv = TRUE,
                                    cv_group_by_site = TRUE){
  # Y: N x P binary response matrix
  # X: N x L environment matrix (first column should be 1s for intercept)
  # Tr: P x K species traits matrix
  # alpha: Elastic net mixing parameter (1 for lasso, 0 for ridge)
  # lambda: "lambda.min"/"lambda.1se" when use_cv = TRUE, or a numeric value when use_cv = FALSE
  # use_cv: if TRUE, choose coefficients from cv.glmnet; if FALSE, fit glmnet at the fixed lambda
  # cv_group_by_site: when use_cv = TRUE, group the CV folds by site so all P
  #   pseudo-likelihood rows of a site stay in one fold (TRUE = default, the
  #   statistically correct choice here). FALSE uses glmnet's random row-level
  #   folds and exists only to reproduce the earlier (leaky) baseline.
  
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
  
  # 4. Fit Penalized Logistic Regression
  # The design matrix glm_X already includes the intercept as its first column 
  # (assuming the first column of X is 1s).
  # We should not penalize the global intercept term.
  penalty_factor = rep(1, n_params)
  penalty_factor[1] = 0 
  
  if (use_cv) {
    if (cv_group_by_site) {
      # Group the cross-validation folds by site. The pseudo-likelihood design has
      # P rows per site (rows are ordered site-by-site: rows (s-1)*P + 1 ... s*P
      # belong to site s) that are dependent by construction -- they share x_s, and
      # each species' presence enters the neighbour terms of the others. Random
      # row-level folds would therefore leak information between train and test and
      # make CV over-optimistic. Putting all P rows of a site in the same fold keeps
      # the folds independent at the site level (sites are generated independently).
      n_folds = min(10, N)
      site_fold = sample(rep(seq_len(n_folds), length.out = N))
      fold_id = rep(site_fold, each = P)

      penalized_fit = cv.glmnet(glm_X, glm_Y, family = "binomial", intercept = FALSE,
                                penalty.factor = penalty_factor, alpha = alpha,
                                foldid = fold_id)
    } else {
      # Default random row-level folds: leaky for these dependent pseudo-likelihood
      # rows, kept only to reproduce the earlier baseline.
      penalized_fit = cv.glmnet(glm_X, glm_Y, family = "binomial", intercept = FALSE,
                                penalty.factor = penalty_factor, alpha = alpha)
    }
  } else {
    if (!is.numeric(lambda)) {
      stop("When use_cv = FALSE, lambda must be a numeric value.")
    }
    penalized_fit = glmnet(glm_X, glm_Y, family = "binomial", intercept = FALSE,
                           penalty.factor = penalty_factor, alpha = alpha,
                           lambda = lambda)
  }
  
  # 5. Extract and reshape parameters
  est_coefs = as.numeric(coef(penalized_fit, s = lambda))
  
  # Note: coef(cv_fit) returns a vector that includes an explicit Intercept term at the start.
  # Since we used intercept = FALSE, this first element is always 0.
  # We drop it to match our glm_X columns.
  est_coefs = est_coefs[-1]
  
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
    penalized_model = penalized_fit,
    cv_model = if (use_cv) penalized_fit else NULL,
    use_cv = use_cv,
    alpha = alpha,
    lambda = lambda
  ))
}