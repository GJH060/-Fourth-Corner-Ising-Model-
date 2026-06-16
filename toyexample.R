rm(list = ls())
library(tidyverse)
library(glmnet)
source("Code/FCIR/generate_without_sparsity.r")
source("Code/FCIR/generate_FCIR.r")
source("Code/FCIR/estimate_penalized_FCIR.r")
source("Code/FCIR/estimate_FCIR.r")


##------------------------
#' # Simulate multivariate binary data from FCIR model
##------------------------
Ns <- 800
Ps <- 60

# Global fixed parameters
L = 3          
K = 2          
seed = 42      


# generate_fully_dense_fcir_data(N = Ns, 
#                                P = Ps, 
#                                L = L, 
#                                K = K, 
#                                B_reps = 1, 
#                                seed = seed, 
#                                filename = "testexample.RData")
generate_dense_fcir_data(N = Ns, 
                         P = Ps, 
                         L = L, 
                         K = K, 
                         B_reps = 1, 
                         seed = seed, 
                         filename = "testexample.RData")
load("testexample.RData")


##------------------------
#' # Fit various forms of the FCIR model 
##------------------------
fit_FCIR_unpen <- estimate_unpenalized_FCIR(Y = Y[,,1],
                                            X = X,
                                            Tr = Tr)


fit_FCIR_ridgelambda <- estimate_penalized_FCIR(Y = Y[,,1],
                                                 X = X,
                                                 Tr = Tr,
                                                 alpha = 0,
                                                 lambda = 0.5 / (Ns * Ps),
                                                 use_cv = FALSE)

fit_FCIR_lassolambda <- estimate_penalized_FCIR(Y = Y[,,1],
                                                X = X,
                                                Tr = Tr,
                                                alpha = 1,
                                                lambda = 0.5 / (Ns * Ps),
                                                use_cv = FALSE)

fit_FCIR_lassolambdacv <- estimate_penalized_FCIR(Y = Y[,,1],
                                                X = X,
                                                Tr = Tr,
                                                alpha = 1,
                                                lambda = "lambda.min",
                                                cv_group_by_site = TRUE,
                                                use_cv = TRUE)


##------------------------
#' # Fit Assess results
##------------------------
data.frame(true = beta_0, 
           unpen = fit_FCIR_unpen$beta_0, 
           ridge = fit_FCIR_ridgelambda$beta_0,
           lasso = fit_FCIR_lassolambda$beta_0,
           lassolambdacv = fit_FCIR_lassolambdacv$beta_0)
data.frame(true = alpha_0, 
           unpen = fit_FCIR_unpen$alpha_0, 
           ridge = fit_FCIR_ridgelambda$alpha_0,
           lasso = fit_FCIR_lassolambda$alpha_0,
           lassolambdacv = fit_FCIR_lassolambdacv$alpha_0)
data.frame(true = B_mat %>% as.vector, 
           unpen = fit_FCIR_unpen$B_mat %>% as.vector,
           ridge = fit_FCIR_ridgelambda$B_mat %>% as.vector,
           lasso = fit_FCIR_lassolambda$B_mat %>% as.vector,
           lassolambdacv = fit_FCIR_lassolambdacv$B_mat %>% as.vector)
data.frame(true = A_mat %>% as.vector, 
           unpen = fit_FCIR_unpen$A_mat %>% as.vector,
           ridge = fit_FCIR_ridgelambda$A_mat %>% as.vector,
           lasso = fit_FCIR_lassolambda$A_mat %>% as.vector,
           lassolambdacv = fit_FCIR_lassolambdacv$A_mat %>% as.vector)


