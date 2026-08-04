rm(list = ls())
library(tidyverse)
library(glmnet)
source("Code/FCIR_M/generate_FCIR_M.r")
source("Code/FCIR_M/estimate_FCIR_M.r")
source("Code/FCIR_M/estimate_adaptive_lasso_FCIR_M.r")


##------------------------
#' # Simulate multivariate binary data from FCIR_M model 
##------------------------
Ns <- 50
Ps <- 20

# Global fixed parameters
L = 3 # Number of covariates
K = 2 # Number of species traits           
seed = 2026


generate_fcir_M_data(N = Ns, 
                     P = Ps, 
                     L = L, 
                     K = K, 
                     B_reps = 1, 
                     seed = seed, 
                     filename = "testexample.RData")
load("testexample.RData")


beta_0
B_mat
Theta_int
t(beta_0 + tcrossprod(B_mat, Tr)) 


##------------------------
#' # Fit various forms of the FCIR model or a special case of it
##------------------------
fit_unpen_nostan <- estimate_unpenalized_FCIR_M(Y = Y[,,1],
                                                X = X,
                                                Tr = Tr,
                                                standardize = FALSE)

fit_adaptivelasso <- estimate_adaptive_lasso_FCIR_M(Y = Y[,,1],
                                                    X = X,
                                                    Tr,
                                                    gamma = 1, 
                                                    lambda = "lambda.min",
                                                    use_cv = TRUE,
                                                    cv_group_by_site = TRUE) 

data.frame(true = beta_0, unpen = fit_unpen_nostan$beta_0, adlasso = fit_adaptivelasso$beta_0)
data.frame(true = B_mat %>% as.vector, unpen = fit_unpen_nostan$B_mat %>% as.vector, adlasso = fit_adaptivelasso$B_mat %>% as.vector)
data.frame(true = Theta_int[lower.tri(Theta_int)], 
           unpen = fit_unpen_nostan$Theta_int[lower.tri(fit_unpen_nostan$Theta_int)], 
           adlasso = fit_adaptivelasso$Theta_int[lower.tri(fit_adaptivelasso$Theta_int)]) %>% 
    t %>% 
    dist

table(Theta_int[lower.tri(Theta_int)] != 0,
      fit_adaptivelasso$Theta_int[lower.tri(fit_adaptivelasso$Theta_int)] != 0)
