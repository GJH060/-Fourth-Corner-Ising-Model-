library(glmnet)

estimate_adaptive_lasso_Ising <- function(Y,
                                          gamma = 1,
                                          init = "unpenalized",
                                          lambda = "lambda.min",
                                          use_cv = TRUE,
                                          symmetrize = "and",
                                          nfolds = 10) {
  # Adaptive lasso on the Ising pseudo-likelihood: one penalized logistic
  # regression per node, with weights from an initial node-wise fit.
  # Thresholds enter as the (unpenalized) glmnet intercept.

  if (!use_cv && !is.numeric(lambda)) {
    stop("When use_cv = FALSE, lambda must be a numeric value.")
  }

  N = nrow(Y)
  P = ncol(Y)
  eps = 1e-6

  init_fit = estimate_nodewise_Ising(Y, init = init, nfolds = nfolds)

  theta_jj = numeric(P)
  Theta_raw = matrix(0, nrow = P, ncol = P)
  sel_lambda = rep(NA_real_, P)
  degenerate = init_fit$degenerate_nodes

  for (j in 1:P) {
    y_j = Y[, j]
    X_j = Y[, -j, drop = FALSE]

    if (degenerate[j]) {
      theta_jj[j] = init_fit$theta_jj[j]
      next
    }

    # Adaptive weights from the initial node-wise slopes.
    beta_init = init_fit$Theta_raw[j, -j]
    beta_init[!is.finite(beta_init)] = 0
    pen_factor = 1 / (abs(beta_init) + eps)^gamma

    fitted_node = tryCatch({
      if (use_cv) {
        foldid = make_stratified_folds(y_j, nfolds)
        if (is.null(foldid)) stop("node cannot be cross-validated")

        fit = cv.glmnet(X_j, y_j, family = "binomial", alpha = 1,
                        penalty.factor = pen_factor, foldid = foldid)
        lam = if (is.character(lambda) &&
                  lambda %in% c("lambda.min", "lambda.1se")) {
          fit[[lambda]]
        } else {
          lambda
        }
      } else {
        fit = glmnet(X_j, y_j, family = "binomial", alpha = 1,
                     penalty.factor = pen_factor, lambda = lambda)
        lam = lambda
      }
      list(cf = as.numeric(coef(fit, s = lam)), lambda = lam)
    }, error = function(e) NULL)

    # A node that cannot be fitted contributes no edges rather than aborting the
    # whole replicate; the count is reported through degenerate_nodes.
    if (is.null(fitted_node)) {
      degenerate[j] = TRUE
      theta_jj[j] = degenerate_threshold(y_j)
      next
    }

    theta_jj[j] = fitted_node$cf[1]
    Theta_raw[j, -j] = fitted_node$cf[-1]
    sel_lambda[j] = fitted_node$lambda
  }

  init_fit$degenerate_nodes = degenerate

  list(
    theta_jj = theta_jj,
    Theta = symmetrize_Theta(Theta_raw, rule = symmetrize),
    Theta_raw = Theta_raw,
    selected_lambda = sel_lambda,
    init_fit = init_fit,
    degenerate_nodes = degenerate,
    gamma = gamma,
    init = init,
    symmetrize = symmetrize,
    use_cv = use_cv
  )
}
