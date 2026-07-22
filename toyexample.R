rm(list = ls())
library(tidyverse)
library(glmnet)
#source("Code/FCIR/generate_without_sparsity.r")
source("Code/FCIR/generate_FCIR.r")
source("Code/FCIR/estimate_FCIR.r")
# source("Code/FCIR/estimate_penalized_FCIR.r")
# source("Code/FCIR/estimate_adaptive_lasso_FCIR.r")

source("Code/FCIR_I/generate_FCIR_I.r")
source("Code/FCIR_I/estimate_FCIR_I.r")
source("Code/FCIR_I/estimate_adaptive_lasso_FCIR_I.r")


##------------------------
#' # Simulate multivariate binary data from FCIR model or a special case of it
##------------------------
Ns <- 1000
Ps <- 30

# Global fixed parameters
L = 3 # Number of covariates
K = 2 # Number of species traits           
seed = 2026


#generate_dense_fcir_data(N = Ns, 
generate_fcir_I_data(N = Ns, 
                     P = Ps, 
                     L = L, 
                     K = K, 
                     B_reps = 1, 
                     seed = seed, 
                     filename = "testexample.RData")
load("testexample.RData")


Beta_mat
A_mat
# t(beta_0 + tcrossprod(B_mat, Tr)) %>% boxplot
# t(beta_0 + tcrossprod(B_mat, Tr)) %>% colMeans
# t(beta_0 + tcrossprod(B_mat, Tr)) %>% apply(., 2, sd)
#apply(beta_0 + B_mat %*% t(Tr), 1, mean)

Delta <- array(0, dim = c(P, P, L))
for(j in 1:P) { for(j_prime in 1:P) {
    Delta[j, j_prime, ] <- alpha_0 + A_mat %*% abs(Tr[j, ] - Tr[j_prime, ])
    } }
apply(Delta, 3, function(x) mean(x[upper.tri(x)]))


##------------------------
#' # Fit various forms of the FCIR model 
##------------------------
# fit_FCIR_unpen_nostan <- estimate_unpenalized_FCIR(Y = Y[,,1],
#                                                    X = X,
#                                                    Tr = Tr,
#                                                    standardize = FALSE)
fit_FCIR_unpen_nostan <- estimate_unpenalized_FCIR_I(Y = Y[,,1],
                                                     X = X,
                                                     Tr = Tr,
                                                     standardize = FALSE)

# fit_FCIR_unpen <- estimate_unpenalized_FCIR(Y = Y[,,1],
#                                             X = X,
#                                             Tr = Tr, 
#                                             standardize = TRUE)
#' The unpenalized fit is on the standardized scale, so we need to rescale the coefficients back to the original scale for comparison with the true parameters. The intercept is adjusted to account for the rescaling of the other coefficients.
# fit_FCIR_unpen$orig_scale_coefficients <- fit_FCIR_unpen$glm_model$coefficients
# fit_FCIR_unpen$orig_scale_coefficients[-1] <- fit_FCIR_unpen$glm_model$coefficients[-1] / fit_FCIR_unpen$getsds  # Rescale by the unpenalized SDs of the Xdes columns
# fit_FCIR_unpen$orig_scale_coefficients[1] <- fit_FCIR_unpen$glm_model$coefficients[1] - sum(fit_FCIR_unpen$glm_model$coefficients[-1] * colMeans(model.matrix(fit_FCIR_unpen_nostan$glm_model)[,-1]) / fit_FCIR_unpen$getsds)  # Adjust intercept for rescaling
# fit_FCIR_unpen$orig_scale_beta0 <- fit_FCIR_unpen$orig_scale_coefficients[1:L]
# fit_FCIR_unpen$orig_scale_B_mat <- matrix(fit_FCIR_unpen$orig_scale_coefficients[(L + 1):(L + L*K)], nrow = L, ncol = K)
# fit_FCIR_unpen$orig_scale_alpha0 <- fit_FCIR_unpen$orig_scale_coefficients[(L + L*K + 1):(2*L + L*K)]
# fit_FCIR_unpen$orig_scale_A_mat <- matrix(fit_FCIR_unpen$orig_scale_coefficients[(2*L + L*K + 1):(2*L + 2*L*K)], nrow = L, ncol = K)
fit_FCIR_unpen <- estimate_unpenalized_FCIR_I(Y = Y[,,1],
                                              X = X,
                                              Tr = Tr,
                                              standardize = TRUE)
#' The unpenalized fit is on the standardized scale, so we need to rescale the coefficients back to the original scale for comparison with the true parameters. 
fit_FCIR_unpen$orig_scale_coefficients <- fit_FCIR_unpen$glm_model$coefficients
fit_FCIR_unpen$orig_scale_coefficients[fit_FCIR_unpen$index_standardized_cols == 1] <- fit_FCIR_unpen$glm_model$coefficients[fit_FCIR_unpen$index_standardized_cols == 1] / fit_FCIR_unpen$getsds[fit_FCIR_unpen$index_standardized_cols == 1]  # Rescale by the unpenalized SDs, for those covariates that were standardized
fit_FCIR_unpen$orig_scale_B_mat <- matrix(fit_FCIR_unpen$orig_scale_coefficients[1:(P*L)], nrow = P, byrow = TRUE)
fit_FCIR_unpen$orig_scale_alpha0 <- fit_FCIR_unpen$orig_scale_coefficients[P*L + 1:(L)]
fit_FCIR_unpen$orig_scale_A_mat <- matrix(fit_FCIR_unpen$orig_scale_coefficients[P*L + L + 1:(L*K)], nrow = L)



# fit_FCIR_lassolambdacv <- estimate_penalized_FCIR(Y = Y[,,1],
#                                                   X = X,
#                                                   Tr = Tr,
#                                                   alpha = 1,
#                                                   lambda = "lambda.min",
#                                                   cv_group_by_site = TRUE,
#                                                   use_cv = TRUE)

# fit_FCIR_adaptivelasso <- estimate_adaptive_lasso_FCIR(Y = Y[,,1],
#                                                        X = X,
#                                                        Tr,
#                                                        gamma = 1, 
#                                                        lambda = "lambda.min",
#                                                        use_cv = TRUE,
#                                                        cv_group_by_site = TRUE) 
fit_FCIR_adaptivelasso <- estimate_adaptive_lasso_FCIR_I(Y = Y[,,1],
                                                         X = X,
                                                         Tr,
                                                         gamma = 1, 
                                                         lambda = "lambda.min",
                                                         use_cv = TRUE,
                                                         cv_group_by_site = TRUE) 


##------------------------
#' # Fit Assess results
##------------------------
beta_0; fit_FCIR_unpen$orig_scale_beta0; fit_FCIR_lassolambdacv$beta_0; fit_FCIR_adaptivelasso$beta_0
B_mat; fit_FCIR_unpen$orig_scale_B; fit_FCIR_lassolambdacv$B_mat; fit_FCIR_adaptivelasso$B_mat

alpha_0; fit_FCIR_unpen$orig_scale_alpha_0; fit_FCIR_lassolambdacv$alpha_0; fit_FCIR_adaptivelasso$alpha_0
A_mat; fit_FCIR_unpen$orig_scale_A; fit_FCIR_lassolambdacv$A_mat; fit_FCIR_adaptivelasso$A_mat




data.frame(true = beta_0, 
           unpen = fit_FCIR_unpen$beta_0, 
           ridge = fit_FCIR_ridgelambda$beta_0,
           lasso = fit_FCIR_lassolambda$beta_0,
           ridgelambdacv = fit_FCIR_ridgelambdacv$beta_0,
           lassolambdacv = fit_FCIR_lassolambdacv$beta_0)
data.frame(true = alpha_0, 
           unpen = fit_FCIR_unpen$alpha_0, 
           ridge = fit_FCIR_ridgelambda$alpha_0,
           lasso = fit_FCIR_lassolambda$alpha_0,
           ridgelambdacv = fit_FCIR_ridgelambdacv$alpha_0,
           lassolambdacv = fit_FCIR_lassolambdacv$alpha_0)
data.frame(true = B_mat %>% as.vector, 
           unpen = fit_FCIR_unpen$B_mat %>% as.vector,
           ridge = fit_FCIR_ridgelambda$B_mat %>% as.vector,
           lasso = fit_FCIR_lassolambda$B_mat %>% as.vector,
           ridgelambdacv = fit_FCIR_ridgelambdacv$B_mat %>% as.vector,
           lassolambdacv = fit_FCIR_lassolambdacv$B_mat %>% as.vector)
data.frame(true = A_mat %>% as.vector, 
           unpen = fit_FCIR_unpen$A_mat %>% as.vector,
           ridge = fit_FCIR_ridgelambda$A_mat %>% as.vector,
           lasso = fit_FCIR_lassolambda$A_mat %>% as.vector,
           ridgelambdacv = fit_FCIR_ridgelambdacv$A_mat %>% as.vector,
           lassolambdacv = fit_FCIR_lassolambdacv$A_mat %>% as.vector)


