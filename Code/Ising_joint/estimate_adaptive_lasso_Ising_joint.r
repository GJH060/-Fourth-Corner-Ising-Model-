library(glmnet)

estimate_adaptive_lasso_Ising_joint <- function(Y,
                                                gamma = 1,
                                                init = "unpenalized",
                                                lambda = "lambda.min",
                                                use_cv = TRUE,
                                                cv_group_by_site = TRUE,
                                                nfolds = 10) {
  # Adaptive lasso on the joint Ising PL: one parameter per undirected edge.
  # Species thresholds (first P columns) are unpenalized.

  if (!init %in% c("unpenalized", "ridge")) {
    stop("init must be 'unpenalized' or 'ridge'.")
  }
  if (!use_cv && !is.numeric(lambda)) {
    stop("When use_cv = FALSE, lambda must be a numeric value.")
  }

  N = nrow(Y)
  P = ncol(Y)
  n_edges = P * (P - 1) / 2
  eps = 1e-6

  make_site_folds <- function() {
    if (!cv_group_by_site) return(NULL)
    n_folds = min(nfolds, N)
    site_fold = sample(rep(seq_len(n_folds), length.out = N))
    rep(site_fold, each = P)
  }

  init_fit = estimate_unpenalized_Ising_joint(Y = Y, standardize = TRUE)
  getsds = init_fit$getsds
  index_std = init_fit$index_standardized_cols
  Xdes = model.matrix(init_fit$glm_model)
  glm_Y = as.numeric(init_fit$glm_model$y)

  if (init == "unpenalized") {
    beta_init = as.numeric(coef(init_fit$glm_model))
  } else {
    foldid = make_site_folds()
    ridge_fit = if (is.null(foldid)) {
      cv.glmnet(Xdes, glm_Y, family = "binomial", alpha = 0,
                intercept = FALSE, standardize = FALSE,
                nfolds = min(nfolds, N))
    } else {
      cv.glmnet(Xdes, glm_Y, family = "binomial", alpha = 0,
                intercept = FALSE, standardize = FALSE, foldid = foldid)
    }
    beta_init = as.numeric(coef(ridge_fit, s = "lambda.min"))
    if (length(beta_init) == ncol(Xdes) + 1) beta_init = beta_init[-1]
  }
  beta_init[!is.finite(beta_init)] = 0

  pen_factor = 1 / (abs(beta_init) + eps)^gamma
  pen_factor[seq_len(P)] = 0

  if (use_cv) {
    ad_fit = cv.glmnet(Xdes, glm_Y, family = "binomial",
                       intercept = FALSE, standardize = FALSE,
                       alpha = 1, penalty.factor = pen_factor,
                       foldid = make_site_folds())
    sel_lambda = if (is.character(lambda) &&
                     lambda %in% c("lambda.min", "lambda.1se")) {
      ad_fit[[lambda]]
    } else {
      lambda
    }
  } else {
    ad_fit = glmnet(Xdes, glm_Y, family = "binomial",
                    intercept = FALSE, standardize = FALSE,
                    alpha = 1, penalty.factor = pen_factor,
                    lambda = lambda)
    sel_lambda = lambda
  }

  est_coefs = as.numeric(coef(ad_fit, s = sel_lambda))
  if (length(est_coefs) == ncol(Xdes) + 1) {
    est_coefs = est_coefs[-1]
  }

  final_coefs = est_coefs
  std_idx = which(index_std == 1)
  final_coefs[std_idx] = est_coefs[std_idx] / getsds[std_idx]

  hat_theta_jj = final_coefs[seq_len(P)]
  hat_theta_int_vec = final_coefs[(P + 1):(P + n_edges)]

  hat_Theta = matrix(0, nrow = P, ncol = P)
  hat_Theta[upper.tri(hat_Theta)] = hat_theta_int_vec
  hat_Theta[lower.tri(hat_Theta)] = t(hat_Theta)[lower.tri(hat_Theta)]

  list(
    theta_jj = hat_theta_jj,
    Theta = hat_Theta,
    selected_lambda = sel_lambda,
    init_fit = init_fit,
    adaptive_model = ad_fit,
    penalty_factor = pen_factor,
    gamma = gamma,
    init = init,
    use_cv = use_cv
  )
}
