estimate_unpenalized_FCIR_M <- function(Y,
                                        X,
                                        Tr,
                                         standardize = FALSE,
                                        returnX_only = FALSE,
                                        sparse = FALSE){
  # Y: N x P binary response matrix
  # X: N x L environment matrix (first column is 1s for intercept)
  # Tr: P x K species traits matrix
  # sparse: build the edge block as a dgCMatrix instead of a dense matrix.
  #   The dense path materialises two N*P x n_edges matrices via outer(), which
  #   is ~15 GB of temporaries at N = 29,983 / P = 30. The edge block is in fact
  #   very sparse -- row (s, j) is non-zero only at edges incident to j, and only
  #   where that neighbour is actually present -- so the sparse path holds
  #   sum_s R_s * (P-1) non-zeros instead of N*P*n_edges cells. Results are
  #   identical; only the storage differs. Default FALSE so the simulation
  #   pipeline is untouched.

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

  if (sparse) {
    if (standardize) {
      stop("sparse = TRUE does not support standardize = TRUE: scale() centers the
 design and destroys sparsity. Scale by SD outside the function instead, which is
 what estimate_adaptive_lasso_FCIR_M() already does.")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("sparse = TRUE requires the Matrix package.")

    # Column e of the edge block corresponds to the pair (edge_i[e], edge_j[e]).
    # It is non-zero on exactly two families of rows:
    #   focal = edge_i[e], carrying the value Y[s, edge_j[e]]
    #   focal = edge_j[e], carrying the value Y[s, edge_i[e]]
    # Iterating over the n_edges columns keeps every inner operation vectorised.
    row_list = vector("list", n_edges)
    col_list = vector("list", n_edges)
    val_list = vector("list", n_edges)
    for (e in seq_len(n_edges)) {
      a = edge_i[e]; b = edge_j[e]
      s_b = which(Y[, b] != 0)          # sites where the neighbour b is present
      s_a = which(Y[, a] != 0)          # sites where the neighbour a is present
      rows = c((s_b - 1L) * P + a, (s_a - 1L) * P + b)
      row_list[[e]] = rows
      col_list[[e]] = rep.int(e, length(rows))
      val_list[[e]] = c(Y[s_b, b], Y[s_a, a])
    }
    comp3_theta_int = Matrix::sparseMatrix(i = unlist(row_list),
                                           j = unlist(col_list),
                                           x = as.numeric(unlist(val_list)),
                                           dims = c(N * P, n_edges))
    glm_X = cbind(Matrix::Matrix(comp1_beta0, sparse = TRUE),
                  Matrix::Matrix(comp2_B, sparse = TRUE),
                  comp3_theta_int)

    # apply(., 2, sd) would densify column by column; compute the same SDs from
    # the sparse column sums instead.
    n_row = nrow(glm_X)
    col_mean = Matrix::colMeans(glm_X)
    col_sqsum = Matrix::colSums(glm_X^2)
    getsds <- sqrt(pmax(0, (col_sqsum - n_row * col_mean^2) / (n_row - 1)))
    names(getsds) <- colnames(glm_X)

    if (returnX_only) {
      return(glm_X)
    }
    # glm() cannot take a sparse design, so densify for the initial fit only and
    # drop the dense copy as soon as it has been consumed.
    glm_X <- as.matrix(glm_X)
  } else {
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

