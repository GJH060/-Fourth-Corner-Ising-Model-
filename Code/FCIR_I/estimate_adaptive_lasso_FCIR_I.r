library(glmnet)

estimate_adaptive_lasso_FCIR_I <- function(Y, X, Tr,
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

  make_site_folds <- function() {
    if (!cv_group_by_site) {
      return(NULL)
    }

    n_folds <- min(10, N)
    site_fold <- sample(rep(seq_len(n_folds), length.out = N))
    rep(site_fold, each = P)
  }

  # Column layout: vec_rows(Beta_mat) | alpha_0 | vec(A).
  # Manual standardization + resolve (mirrors estimate_adaptive_lasso_FCIR.r).
  # FCIR_I has species-specific intercepts inside the design and is fit with
  # intercept = FALSE, so columns are scaled by their SD but NOT centered, to
  # keep the linear predictor intact without a single global intercept.
  init_fit <- estimate_unpenalized_FCIR_I(Y = Y, 
                                          X = X, 
                                          Tr = Tr, 
                                          standardize = TRUE)
  # getsds <- init_fit$getsds
  # Xdes_raw <- estimate_unpenalized_FCIR_I(Y = Y, X = X, Tr = Tr, returnX_only = TRUE)
  # index_intercept_cols <- seq(1, P * L, by = L)
  # Xdes <- Xdes_raw
  # Xdes[,-index_intercept_cols] <- sweep(Xdes_raw[,-index_intercept_cols], 2, getsds[-index_intercept_cols], "/")   # standardized design (SD scaling only)
  glm_Y <- as.numeric(init_fit$glm_model$y)

  # Initial estimator on the standardized scale for the adaptive weights:
  # rescaling column c by 1/sd_c multiplies its unpenalized coefficient by sd_c.
  # beta_init <- as.numeric(coef(init_fit$glm_model)) 
  # beta_init[-index_intercept_cols] <- beta_init[-index_intercept_cols] * getsds[-index_intercept_cols]   # rescale to raw scale
  beta_init <- as.numeric(coef(init_fit$glm_model))
  beta_init[!is.finite(beta_init)] <- 0
  Xdes <- model.matrix(init_fit$glm_model)
  
  eps <- 1e-6
  pen_factor <- 1 / (abs(beta_init) + eps)^gamma
  if (length(custom_penalty_factor) != 1 &&
      length(custom_penalty_factor) != length(pen_factor)) {
    stop(sprintf("custom_penalty_factor must have length 1 or %d, got %d.",
                 length(pen_factor), length(custom_penalty_factor)))
  }
  pen_factor <- pen_factor * custom_penalty_factor

  # Separate out the species-specific intercepts: the first column of each
  # species block (the x_s[1] = 1 term) is left unpenalized, mirroring how FCIR
  # leaves the global intercept beta_0[1] unpenalized.
  intercept_idx <- seq(1, P * L, by = L)
  pen_factor[intercept_idx] <- 0

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
  final_coefs <- est_coefs
  final_coefs[init_fit$index_standardized_cols == 1] <- est_coefs[init_fit$index_standardized_cols == 1] / init_fit$getsds[init_fit$index_standardized_cols == 1]
  
  
  idx <- 1
  hat_Beta_vec <- final_coefs[idx:(idx + P * L - 1)]; idx <- idx + P * L
  hat_alpha_0  <- final_coefs[idx:(idx + L - 1)]; idx <- idx + L
  hat_A_vec    <- final_coefs[idx:(idx + L * K - 1)]

  list(
    Beta_mat = t(matrix(hat_Beta_vec, nrow = L, ncol = P)),
    alpha_0 = hat_alpha_0,
    A_mat = matrix(hat_A_vec, nrow = L, ncol = K),
    adaptive_model = ad_fit,
    cv_model = if (use_cv) ad_fit else NULL,
    init_model = init_fit,
    penalty_factor = pen_factor,
    selected_lambda = sel_lambda,
    gamma = gamma,
    use_cv = use_cv
  )
}
