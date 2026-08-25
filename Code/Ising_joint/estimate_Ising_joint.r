library(glmnet)

# Joint pseudo-likelihood for a plain Ising model with one parameter per
# undirected edge (upper-triangle parameterization), analogous to FCIR_M's
# Theta_int block.
#
#   logit P(y_j = 1 | y_-j) = theta_jj + sum_{j' != j} theta_jj' y_j'
#
# Stacking over sites and species yields an NP-row logistic regression with
# columns [species thresholds (P) | edge indicators (P*(P-1)/2)].

build_ising_pl_design_joint <- function(Y) {
  N = nrow(Y)
  P = ncol(Y)
  n_edges = P * (P - 1) / 2

  pair_idx_mat = matrix(0, nrow = P, ncol = P)
  pair_idx_mat[upper.tri(pair_idx_mat)] = 1:n_edges
  pair_idx_mat[lower.tri(pair_idx_mat)] = t(pair_idx_mat)[lower.tri(pair_idx_mat)]

  n_obs = N * P
  n_params = P + n_edges
  glm_Y = numeric(n_obs)
  glm_X = matrix(0, nrow = n_obs, ncol = n_params)

  row_idx = 1
  for (s in 1:N) {
    y_s = Y[s, ]
    for (j in 1:P) {
      glm_Y[row_idx] = y_s[j]
      glm_X[row_idx, j] = 1

      present = which(y_s == 1)
      present = present[present != j]
      if (length(present) > 0) {
        glm_X[row_idx, P + pair_idx_mat[j, present]] = 1
      }
      row_idx = row_idx + 1
    }
  }

  list(glm_Y = glm_Y, glm_X = glm_X, N = N, P = P, n_edges = n_edges,
       pair_idx_mat = pair_idx_mat)
}

estimate_unpenalized_Ising_joint <- function(Y,
                                             standardize = FALSE,
                                             returnX_only = FALSE) {
  des = build_ising_pl_design_joint(Y)
  glm_Y = des$glm_Y
  glm_X = des$glm_X
  P = des$P
  n_edges = des$n_edges

  getsds = apply(glm_X, 2, sd)
  index_standardized_cols = numeric(ncol(glm_X))

  if (standardize) {
    # Leave species-threshold columns unscaled; scale edge columns by SD only
    # (no centering), mirroring FCIR_I / FCIR_M.
    edge_idx = (P + 1):(P + n_edges)
    valid_idx = edge_idx[is.finite(getsds[edge_idx]) & getsds[edge_idx] > 0]
    scale_factors = rep(1, ncol(glm_X))
    scale_factors[valid_idx] = getsds[valid_idx]
    glm_X = sweep(glm_X, 2, scale_factors, "/")
    index_standardized_cols[valid_idx] = 1
  }

  if (returnX_only) {
    return(glm_X)
  }

  logistic_reg = glm(glm_Y ~ glm_X + 0, family = binomial,
                     control = glm.control(maxit = 100))
  est_coefs = logistic_reg$coefficients
  est_coefs[!is.finite(est_coefs)] = 0

  hat_theta_jj = est_coefs[1:P]
  hat_theta_int_vec = est_coefs[(P + 1):(P + n_edges)]

  hat_Theta = matrix(0, nrow = P, ncol = P)
  hat_Theta[upper.tri(hat_Theta)] = hat_theta_int_vec
  hat_Theta[lower.tri(hat_Theta)] = t(hat_Theta)[lower.tri(hat_Theta)]

  list(
    theta_jj = hat_theta_jj,
    Theta = hat_Theta,
    standardize = standardize,
    getsds = getsds,
    index_standardized_cols = index_standardized_cols,
    glm_model = logistic_reg,
    design = des
  )
}
