library(glmnet)

estimate_adaptive_lasso_FCIR_M <- function(Y, X, Tr,
                                           gamma = 1,
                                           init = "unpenalized",
                                           lambda = "lambda.min",
                                           use_cv = TRUE,
                                           cv_group_by_site = TRUE,
                                           custom_penalty_factor = 1) {
  if (init != "unpenalized") {
    stop("Only init = 'unpenalized' is currently supported.")
  }

  if (!use_cv && !is.numeric(lambda)) {
    stop("When use_cv = FALSE, lambda must be a numeric value.")
  }

  N <- nrow(Y)
  P <- ncol(Y)
  L <- ncol(X)
  K <- ncol(Tr)
  n_edges <- P * (P - 1) / 2

  make_site_folds <- function() {
    if (!cv_group_by_site) {
      return(NULL)
    }

    n_folds <- min(10, N)
    site_fold <- sample(rep(seq_len(n_folds), length.out = N))
    rep(site_fold, each = P)
  }

  # Column layout: beta_0 (L) | vec(B) (L*K) | theta_int (n_edges).
  # Manual standardization + resolve (mirrors estimate_adaptive_lasso_FCIR_I.r).
  # FCIR_M shares FCIR's single global intercept beta_0[1] (design column 1 is a
  # constant ones column). Columns are scaled by SD but NOT centered, the
  # constant intercept column is left unscaled (SD = 0), and it is fit with
  # intercept = FALSE so beta_0[1] is recovered as an unpenalized coefficient.
  init_fit <- estimate_unpenalized_FCIR_M(Y = Y, X = X, Tr = Tr, standardize = FALSE)
  getsds <- init_fit$getsds
  scale_sds <- getsds
  scale_sds[scale_sds == 0] <- 1   # do not scale constant columns (e.g. intercept)
  Xdes_raw <- estimate_unpenalized_FCIR_M(Y = Y, X = X, Tr = Tr, returnX_only = TRUE)
  Xdes <- sweep(Xdes_raw, 2, scale_sds, "/")   # standardized design (SD scaling only)
  glm_Y <- as.numeric(init_fit$glm_model$y)

  # Initial estimator on the standardized scale for the adaptive weights:
  # rescaling column c by 1/sd_c multiplies its unpenalized coefficient by sd_c.
  beta_init <- as.numeric(coef(init_fit$glm_model)) * scale_sds
  beta_init[!is.finite(beta_init)] <- 0

  eps <- 1e-6
  pen_factor <- 1 / (abs(beta_init) + eps)^gamma
  if (length(custom_penalty_factor) != 1 &&
      length(custom_penalty_factor) != length(pen_factor)) {
    stop(sprintf("custom_penalty_factor must have length 1 or %d, got %d.",
                 length(pen_factor), length(custom_penalty_factor)))
  }
  pen_factor <- pen_factor * custom_penalty_factor

  # Separate out the global intercept: beta_0[1] (design column 1, the x_s[1] = 1
  # term) is left unpenalized, mirroring how FCIR leaves beta_0[1] unpenalized.
  pen_factor[1] <- 0

  # Adaptive lasso on the manually standardized design (standardize = FALSE).
  if (use_cv) {
    ad_fit <- cv.glmnet(Xdes,
                        glm_Y,
                        family = "binomial",
                        intercept = FALSE,
                        standardize = FALSE,
                        alpha = 1,
                        penalty.factor = pen_factor,
                        foldid = make_site_folds())

    sel_lambda <- if (is.character(lambda) && lambda %in% c("lambda.min", "lambda.1se")) {
      ad_fit[[lambda]]
    } else {
      lambda
    }
  } else {
    ad_fit <- glmnet(Xdes,
                     glm_Y,
                     family = "binomial",
                     intercept = FALSE,
                     standardize = FALSE,
                     alpha = 1,
                     penalty.factor = pen_factor,
                     lambda = lambda)
    sel_lambda <- lambda
  }

  est_coefs <- as.numeric(coef(ad_fit, s = sel_lambda))
  if (length(est_coefs) == ncol(Xdes) + 1) {
    est_coefs <- est_coefs[-1]   # drop the (zero) glmnet intercept slot
  }

  # Resolve: standardized coefficient / sd = coefficient on the raw scale.
  final_coefs <- est_coefs / scale_sds

  idx <- 1
  hat_beta_0        <- final_coefs[idx:(idx + L - 1)]; idx <- idx + L
  hat_B_vec         <- final_coefs[idx:(idx + L * K - 1)]; idx <- idx + L * K
  hat_theta_int_vec <- final_coefs[idx:(idx + n_edges - 1)]

  hat_Theta_int <- matrix(0, nrow = P, ncol = P)
  hat_Theta_int[upper.tri(hat_Theta_int)] <- hat_theta_int_vec
  hat_Theta_int[lower.tri(hat_Theta_int)] <- t(hat_Theta_int)[lower.tri(hat_Theta_int)]

  list(
    beta_0 = hat_beta_0,
    B_mat = matrix(hat_B_vec, nrow = L, ncol = K),
    Theta_int = hat_Theta_int,
    adaptive_model = ad_fit,
    cv_model = if (use_cv) ad_fit else NULL,
    init_model = init_fit,
    penalty_factor = pen_factor,
    selected_lambda = sel_lambda,
    gamma = gamma,
    use_cv = use_cv
  )
}
