library(glmnet)
estimate_adaptive_lasso_FCIR <- function(Y, X, Tr,
                                         gamma = 1,
                                         init = c("ridge", "unpenalized"),
                                         lambda = "lambda.min",
                                         use_cv = TRUE,
                                         cv_group_by_site = TRUE,
                                         init_lambda = "lambda.min") {
  init <- match.arg(init)
  N <- nrow(Y); P <- ncol(Y); L <- ncol(X); K <- ncol(Tr)

  # 1. Pre-compute trait differences
  Delta <- array(0, dim = c(P, P, K))
  for (j in 1:P) for (j_prime in 1:P) Delta[j, j_prime, ] <- abs(Tr[j, ] - Tr[j_prime, ])

  # 2. Build the pseudo-likelihood design (rows ordered site-by-site)
  n_obs <- N * P
  n_params <- 2 * L + 2 * L * K
  glm_Y <- numeric(n_obs)
  glm_X <- matrix(0, nrow = n_obs, ncol = n_params)
  row_idx <- 1
  for (s in 1:N) {
    x_s <- X[s, ]; y_s <- Y[s, ]; R_s <- sum(y_s)
    for (j in 1:P) {
      glm_Y[row_idx] <- y_s[j]
      t_j <- Tr[j, ]
      comp2_B <- kronecker(t_j, x_s)
      neighbor_sum <- R_s - y_s[j]
      comp3_alpha0 <- neighbor_sum * x_s
      w_sj <- numeric(K)
      for (j_prime in 1:P) if (j_prime != j && y_s[j_prime] == 1) w_sj <- w_sj + Delta[j, j_prime, ]
      comp4_A <- kronecker(w_sj, x_s)
      glm_X[row_idx, ] <- c(x_s, comp2_B, comp3_alpha0, comp4_A)
      row_idx <- row_idx + 1
    }
  }
  # Drop the all-ones intercept column; glmnet fits its own intercept.
  Xdes <- glm_X[, -1, drop = FALSE]

  # Xdes column layout (length n_params - 1): beta_0[2:L] | vec(B) | alpha_0 | vec(A).
  # The final coef() vector restores beta_0[1] as the glmnet intercept (see step 6).

  # Site-grouped folds: all P rows of a site stay in one fold.
  make_site_folds <- function() {
    n_folds <- min(10, N)
    site_fold <- sample(rep(seq_len(n_folds), length.out = N))
    rep(site_fold, each = P)
  }

  # 3. Initial estimator beta_init (excludes the glmnet/glm intercept)
  if (init == "ridge") {
    if (use_cv && cv_group_by_site) {
      init_fit <- cv.glmnet(Xdes, glm_Y, family = "binomial", intercept = TRUE,
                            alpha = 0, foldid = make_site_folds())
    } else {
      init_fit <- cv.glmnet(Xdes, glm_Y, family = "binomial", intercept = TRUE, alpha = 0)
    }
    beta_init <- as.numeric(coef(init_fit, s = init_lambda))[-1]
  } else {
    init_fit <- glm(glm_Y ~ Xdes + 1, family = binomial, control = glm.control(maxit = 100))
    beta_init <- as.numeric(coef(init_fit))[-1]
  }
  # Aliased/NA initial coefficients (rank deficiency) -> treat as ~0 (heavy penalty).
  beta_init[is.na(beta_init)] <- 0

  # 4. Adaptive penalty on ALL non-intercept coefficients .
  eps <- 1e-6
  pen_factor <- 1 / (abs(beta_init) + eps)^gamma
  pen_factor <- pen_factor / mean(pen_factor)
  pen_factor <- pmin(pmax(pen_factor, 1e-3), 1e3)

  # 5. Adaptive lasso fit (alpha = 1)
  if (use_cv) {
    if (cv_group_by_site) {
      ad_fit <- cv.glmnet(Xdes, glm_Y, family = "binomial", intercept = TRUE,
                          alpha = 1, penalty.factor = pen_factor, foldid = make_site_folds())
    } else {
      ad_fit <- cv.glmnet(Xdes, glm_Y, family = "binomial", intercept = TRUE,
                          alpha = 1, penalty.factor = pen_factor)
    }
    sel_lambda <- if (lambda %in% c("lambda.min", "lambda.1se")) ad_fit[[lambda]] else lambda
  } else {
    if (!is.numeric(lambda)) stop("When use_cv = FALSE, lambda must be a numeric value.")
    ad_fit <- glmnet(Xdes, glm_Y, family = "binomial", intercept = TRUE,
                     alpha = 1, penalty.factor = pen_factor, lambda = lambda)
    sel_lambda <- lambda
  }

  # 6. Extract coefficients. coef() returns [intercept, Xdes coefs]; with the
  # ones-column dropped, that intercept IS beta_0[1], so the full vector maps
  # directly onto (beta_0, vec(B), alpha_0, vec(A)).
  est_coefs <- as.numeric(coef(ad_fit, s = lambda))

  idx <- 1
  hat_beta_0  <- est_coefs[idx:(idx + L - 1)]; idx <- idx + L
  hat_B_vec   <- est_coefs[idx:(idx + L * K - 1)]; idx <- idx + L * K
  hat_alpha_0 <- est_coefs[idx:(idx + L - 1)]; idx <- idx + L
  hat_A_vec   <- est_coefs[idx:(idx + L * K - 1)]

  list(
    beta_0 = hat_beta_0,
    B_mat = matrix(hat_B_vec, nrow = L, ncol = K),
    alpha_0 = hat_alpha_0,
    A_mat = matrix(hat_A_vec, nrow = L, ncol = K),
    adaptive_model = ad_fit,
    cv_model = if (use_cv) ad_fit else NULL,
    init_model = init_fit,
    beta_init = beta_init,
    penalty_factor = pen_factor,
    selected_lambda = sel_lambda,
    gamma = gamma,
    init = init,
    use_cv = use_cv
  )
}
