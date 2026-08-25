estimate_unpenalized_FCIR_M <- function(Y, 
                                        X, 
                                        Tr,
                                         standardize = FALSE,
                                        returnX_only = FALSE){
  # Y: N x P binary response matrix
  # X: N x L environment matrix (first column is 1s for intercept)
  # Tr: P x K species traits matrix
  
  N = nrow(Y)
  P = ncol(Y)
  L = ncol(X)
  K = ncol(Tr)
  n_edges = P * (P - 1) / 2
  
  # Build the stacked pseudo-likelihood design without repeatedly assigning
  # matrix subsets. This avoids an R 4.6 JIT/foreach locked-binding error.
  edge_pairs = which(upper.tri(matrix(FALSE, P, P)), arr.ind = TRUE)
  edge_i = edge_pairs[, 1]
  edge_j = edge_pairs[, 2]

  site_idx = rep(seq_len(N), each = P)
  focal_idx = rep(seq_len(P), times = N)

  glm_Y = as.numeric(t(Y))
  comp1_beta0 = X[site_idx, , drop = FALSE]

  trait_rows = Tr[focal_idx, , drop = FALSE]
  comp2_B = trait_rows[, rep(seq_len(K), each = L), drop = FALSE] *
    comp1_beta0[, rep(seq_len(L), times = K), drop = FALSE]

  comp3_theta_int =
    outer(focal_idx, edge_i, "==") * Y[site_idx, edge_j, drop = FALSE] +
    outer(focal_idx, edge_j, "==") * Y[site_idx, edge_i, drop = FALSE]

  glm_X = cbind(comp1_beta0, comp2_B, comp3_theta_int)
  
  getsds <- apply(glm_X, 2, sd)
  if (standardize) {
    glm_X <- scale(glm_X)
  }

  if (returnX_only) {
    return(glm_X)
  }

  # 3. Fit Unpenalized Logistic Regression
  logistic_reg = glm(glm_Y ~ glm_X + 0, family = binomial)
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
  
  return(list(beta_0 = hat_beta_0,
              B_mat = hat_B_mat,
              Theta_int = hat_Theta_int,
              standardize = standardize,
              getsds = getsds,
              glm_model = logistic_reg))
  }

