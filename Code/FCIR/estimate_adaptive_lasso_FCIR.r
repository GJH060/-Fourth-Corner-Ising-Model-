library(glmnet)
estimate_adaptive_lasso_FCIR <- function(Y, X, Tr,
                                         gamma = 1,
                                         init = "unpenalized",
                                         lambda = "lambda.min",
                                         use_cv = TRUE,
                                         cv_group_by_site = TRUE,
                                         custom_penalty_factor = 1) {
    if (init != "unpenalized") {
        stop("Only init = 'unpenalized' is currently supported.")
        }
    
    if(!use_cv && !is.numeric(lambda)) {
        stop("When use_cv = FALSE, lambda must be a numeric value.")
        }
     
    N <- nrow(Y)
    P <- ncol(Y)
    L <- ncol(X)
    K <- ncol(Tr)
    
    # Site-grouped folds: all P rows of a site stay in one fold.
    make_site_folds <- function() { 
        if(!cv_group_by_site) {
            return(NULL)
            } 
        
        n_folds <- min(10, N)
        site_fold <- sample(rep(seq_len(n_folds), length.out = N))
        return(rep(site_fold, each = P))
        }
    
    # Initial estimator beta_init (excludes the glmnet/glm intercept)
    # Xdes column layout (length n_params - 1): beta_0[2:L] | vec(B) | alpha_0 | vec(A).
    # The final coef() vector restores beta_0[1] as the glmnet intercept (see step 6).
    init_fit <- estimate_unpenalized_FCIR(Y = Y,
                                          X = X,
                                          Tr = Tr,
                                          standardize = TRUE)
    beta_init <- as.numeric(coef(init_fit$glm_model))[-1]
    beta_init[!is.finite(beta_init)] <- 0
    Xdes <- model.matrix(init_fit$glm_model)[, -1, drop = FALSE] 
    glm_Y <- as.vector(init_fit$glm_model$y)
    
    
    # Adaptive-lasso penalty weights for all non-intercept coefficients, optionally rescaled per-coefficient by custom_penalty_factor
    # (layout: beta_0[2:L] | vec(B) | alpha_0 | vec(A); default 1 = no change).
    eps <- 1e-6
    pen_factor <- 1 / (abs(beta_init) + eps)^gamma
    if (length(custom_penalty_factor) != 1 &&
        length(custom_penalty_factor) != length(pen_factor)) {
        stop(sprintf("custom_penalty_factor must have length 1 or %d, got %d.",
                     length(pen_factor), length(custom_penalty_factor)))
    }
    pen_factor <- pen_factor * custom_penalty_factor
    
    
    # Adaptive lasso fit (alpha = 1)
    if (use_cv) {
        ad_fit <- cv.glmnet(Xdes, 
                            glm_Y, 
                            family = "binomial",
                            intercept = TRUE, 
                            standardize = FALSE,
                            alpha = 1, 
                            penalty.factor = pen_factor, 
                            foldid = make_site_folds())
        
        sel_lambda <- if (lambda %in% c("lambda.min", "lambda.1se")) ad_fit[[lambda]] else lambda
        } else {
            ad_fit <- glmnet(Xdes, 
                             glm_Y, 
                             family = "binomial", 
                             intercept = TRUE, 
                             standardize = TRUE,
                             alpha = 1, 
                             penalty.factor = pen_factor, 
                             lambda = lambda)
            sel_lambda <- lambda
        }
    
    # Extract coefficients. coef() returns [intercept, Xdes coefs]; with the ones-column dropped, that intercept IS beta_0[1], so the full vector maps directly onto (beta_0, vec(B), alpha_0, vec(A)).
    est_coefs <- as.numeric(coef(ad_fit, s = lambda))
    final_coefs <- est_coefs
    final_coefs[-1] <- est_coefs[-1] / init_fit$getsds # Rescale slopes by the unpenalized SDs of the Xdes columns   
    getunstandardizedX <-  estimate_unpenalized_FCIR(Y = Y,
                                                     X = X,
                                                     Tr = Tr,
                                                     standardize = FALSE,
                                                     returnX_only = TRUE)
    final_coefs[1] <- est_coefs[1] - sum(est_coefs[-1] * colMeans(getunstandardizedX[,-1,drop=FALSE]) / init_fit$getsds)  # Adjust intercept for rescaling
    rm(est_coefs)
    
    idx <- 1
    hat_beta_0  <- final_coefs[idx:(idx + L - 1)]; idx <- idx + L
    hat_B_vec   <- final_coefs[idx:(idx + L * K - 1)]; idx <- idx + L * K
    hat_alpha_0 <- final_coefs[idx:(idx + L - 1)]; idx <- idx + L
    hat_A_vec   <- final_coefs[idx:(idx + L * K - 1)]
    
    list(
        beta_0 = hat_beta_0,
        B_mat = matrix(hat_B_vec, nrow = L, ncol = K),
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

