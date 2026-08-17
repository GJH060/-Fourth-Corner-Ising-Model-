library(glmnet)

# Node-wise (pseudo-likelihood) estimation of a plain Ising model.
#
# With responses coded {0, 1}, the Ising conditional distribution is
#   logit P(y_j = 1 | y_-j) = theta_jj + sum_{i != j} theta_jj' y_i
# so a logistic regression of column j on the remaining columns recovers the
# threshold as its intercept and the row Theta[j, -j] as its slopes. Fitting
# every node gives an asymmetric Theta that a symmetrization rule then resolves.

symmetrize_Theta <- function(Theta_raw, rule = "and") {
  Theta_sym = (Theta_raw + t(Theta_raw)) / 2

  if (rule == "and") {
    # Keep an edge only when both node-wise fits selected it.
    drop = (Theta_raw == 0) | (t(Theta_raw) == 0)
    Theta_sym[drop] = 0
  } else if (rule != "or") {
    stop("rule must be 'and' or 'or'.")
  }

  diag(Theta_sym) = 0
  Theta_sym
}

# glmnet refuses a binomial response with fewer than 2 observations in a class,
# and cv.glmnet fits on training folds, so that has to hold inside every fold
# too. Random 10-fold assignment can strand the rare class even when the node
# as a whole looks fine, which is what makes the failure intermittent. Folds are
# therefore stratified, and their number is capped by the rare-class count so
# each training set keeps at least 2 of it.
make_stratified_folds <- function(y, nfolds = 10) {
  n_min = min(sum(y == 0), sum(y == 1))

  # k = min(nfolds, n_min) holds out at most ceiling(n_min / k) of the rare
  # class, leaving n_min - ceiling(n_min / k) >= 2 exactly when n_min >= 3.
  if (n_min < 3) return(NULL)
  k = min(nfolds, n_min)

  foldid = integer(length(y))
  for (cls in c(0, 1)) {
    idx = which(y == cls)
    foldid[idx] = sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

# Threshold for a node the neighbourhood regression cannot be fitted on.
degenerate_threshold <- function(y) {
  N = length(y)
  p_hat = min(max(mean(y), 1 / (2 * N)), 1 - 1 / (2 * N))
  log(p_hat / (1 - p_hat))
}

estimate_nodewise_Ising <- function(Y, init = "ridge", nfolds = 10) {
  # init = "unpenalized": plain glm per node. Fast, but the MLE does not exist
  #   under (quasi-)separation, which is common when N is small relative to P
  #   or the response is saturated; coefficients then diverge.
  # init = "ridge": cv.glmnet(alpha = 0) per node. Always finite, so it is the
  #   safer source of adaptive-lasso weights (Zou 2006, Sec. 3.5).

  if (!init %in% c("unpenalized", "ridge")) {
    stop("init must be 'unpenalized' or 'ridge'.")
  }

  N = nrow(Y)
  P = ncol(Y)

  theta_jj = numeric(P)
  Theta_raw = matrix(0, nrow = P, ncol = P)
  degenerate = logical(P)
  separated = logical(P)

  for (j in 1:P) {
    y_j = Y[, j]
    X_j = Y[, -j, drop = FALSE]

    foldid = make_stratified_folds(y_j, nfolds)

    # Too few observations in one class for a fittable neighbourhood: pin the
    # threshold to a clipped log-odds and leave the node's edges at zero.
    if (is.null(foldid)) {
      degenerate[j] = TRUE
      theta_jj[j] = degenerate_threshold(y_j)
      next
    }

    cf = tryCatch({
      if (init == "unpenalized") {
        fit = suppressWarnings(
          glm(y_j ~ X_j, family = binomial, control = glm.control(maxit = 100))
        )
        # glm reports converged = TRUE once the deviance stops moving, which
        # also happens when coefficients run off under separation. Flag it on
        # the fitted probabilities instead.
        pr = fitted(fit)
        separated[j] = any(pr < 1e-8 | pr > 1 - 1e-8)
        out = coef(fit)
        out[!is.finite(out)] = 0
        out
      } else {
        fit = cv.glmnet(X_j, y_j, family = "binomial", alpha = 0,
                        foldid = foldid)
        as.numeric(coef(fit, s = "lambda.min"))
      }
    }, error = function(e) NULL)

    if (is.null(cf)) {
      degenerate[j] = TRUE
      theta_jj[j] = degenerate_threshold(y_j)
      next
    }

    theta_jj[j] = cf[1]
    Theta_raw[j, -j] = cf[-1]
  }

  list(
    theta_jj = theta_jj,
    Theta_raw = Theta_raw,
    Theta = symmetrize_Theta(Theta_raw, rule = "or"),
    init = init,
    degenerate_nodes = degenerate,
    separated_nodes = separated
  )
}
