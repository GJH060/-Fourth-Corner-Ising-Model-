rm(list = ls())
library(tidyverse)
library(glmnet)
source("Code/FCIR_I/generate_FCIR_I.r")
source("Code/FCIR_I/estimate_FCIR_I.r")
source("Code/FCIR_I/estimate_adaptive_lasso_FCIR_I.r")


##------------------------
#' # Simulate multivariate binary data from FCIR model or a special case of it
##------------------------
Ns <- 200
Ps <- 60

# Global fixed parameters
L = 3 # Number of covariates
K = 2 # Number of species traits           
seed = 2026

generate_fcir_I_data(N = Ns, 
                     P = Ps, 
                     L = L, 
                     K = K, 
                     B_reps = 1, 
                     seed = seed, 
                     filename = "testexample.RData")
load("testexample.RData")


Beta_mat
alpha_0
A_mat
colMeans(Y[,,1])
rowMeans(Y[,,1]) %>% table


##------------------------
#' # Fit various forms of the FCIR_I model 
##------------------------
fit_unpen_nostan <- estimate_unpenalized_FCIR_I(Y = Y[,,1],
                                                X = X,
                                                Tr = Tr,
                                                standardize = FALSE)

data.frame(true = alpha_0, unpen = fit_unpen_nostan$alpha_0)
data.frame(true = A_mat %>% as.vector, unpen = fit_unpen_nostan$A_mat %>% as.vector)
data.frame(true = Beta_mat[,-1] %>% as.vector, unpen = fit_unpen_nostan$Beta_mat[,-1] %>% as.vector) 
data.frame(true = Beta_mat[,-1] %>% as.vector, unpen = fit_unpen_nostan$Beta_mat[,-1] %>% as.vector) %>% 
    plot
abline(0,1)


#' fit_unpen <- estimate_unpenalized_FCIR_I(Y = Y[,,1],
#'                                               X = X,
#'                                               Tr = Tr,
#'                                               standardize = TRUE)
#' #' The unpenalized fit is on the standardized scale, so we need to rescale the coefficients back to the original scale for comparison with the true parameters.
#' fit_unpen$orig_scale_coefficients <- fit_unpen$glm_model$coefficients
#' fit_unpen$orig_scale_coefficients[fit_unpen$index_standardized_cols == 1] <- fit_unpen$glm_model$coefficients[fit_unpen$index_standardized_cols == 1] / fit_unpen$getsds[fit_unpen$index_standardized_cols == 1]  # Rescale by the unpenalized SDs, for those covariates that were standardized
#' fit_unpen$orig_scale_B_mat <- matrix(fit_unpen$orig_scale_coefficients[1:(P*L)], nrow = P, byrow = TRUE)
#' fit_unpen$orig_scale_alpha0 <- fit_unpen$orig_scale_coefficients[P*L + 1:(L)]
#' fit_unpen$orig_scale_A_mat <- matrix(fit_unpen$orig_scale_coefficients[P*L + L + 1:(L*K)], nrow = L)


fit_adaptivelasso <- estimate_adaptive_lasso_FCIR_I(Y = Y[,,1],
                                                    X = X,
                                                    Tr,
                                                    gamma = 1, 
                                                    lambda = "lambda.min",
                                                    use_cv = TRUE,
                                                    cv_group_by_site = TRUE) 

data.frame(true = alpha_0, 
           unpen = fit_unpen_nostan$alpha_0, 
           adaptivelasso = fit_adaptivelasso$alpha_0)

data.frame(true = A_mat %>% as.vector,
           unpen = fit_unpen_nostan$A_mat %>% as.vector,
           adaptivelasso = fit_adaptivelasso$A_mat %>% as.vector)

data.frame(true = Beta_mat[,-1] %>% as.vector,
           unpen = fit_unpen_nostan$Beta_mat[,-1] %>% as.vector,
           adaptivelasso = fit_adaptivelasso$Beta_mat[,-1] %>% as.vector)

table(true = as.vector(Beta_mat[,-1]) != 0, 
      adaptivelasso = as.vector(fit_adaptivelasso$Beta_mat[,-1]) != 0)
