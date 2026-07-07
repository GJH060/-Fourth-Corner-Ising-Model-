rm(list = ls())
library(tidyverse)
library(glmnet)
library(glm2)
source("Code/FCIR/generate_without_sparsity.r")
source("Code/FCIR/generate_FCIR.r")
source("Code/FCIR/estimate_penalized_FCIR.r")
source("Code/FCIR/estimate_adaptive_lasso_FCIR.r")
source("Code/FCIR/estimate_FCIR.r")


##------------------------
#' # Simulate multivariate binary data from FCIR model
##------------------------
Ns <- 800
Ps <- 30

# Global fixed parameters
L = 3          
K = 2          
seed = 2026


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


beta_0 + B_mat %*% t(Tr)
apply(beta_0 + B_mat %*% t(Tr), 1, mean)

Delta <- array(0, dim = c(P, P, L))
for(j in 1:P) { for(j_prime in 1:P) {
    Delta[j, j_prime, ] <- alpha_0 + A_mat %*% abs(Tr[j, ] - Tr[j_prime, ])
    } }
apply(Delta, 3, function(x) mean(x[upper.tri(x)]))


##------------------------
#' # Fit various forms of the FCIR model 
##------------------------
fit_FCIR_unpen <- estimate_unpenalized_FCIR(Y = Y[,,1],
                                            X = X,
                                            Tr = Tr)

beta_0; fit_FCIR_unpen$beta_0
B_mat; fit_FCIR_unpen$B_mat
alpha_0; fit_FCIR_unpen$alpha_0
A_mat; fit_FCIR_unpen$A_mat


MM <- fit_FCIR_unpen$glm_model %>% model.matrix
summary(abs(MM[,1:(L + L*K)] %*% c(beta_0, B_mat %>% as.vector)))
summary(abs(MM[,(L + L*K + 1):(2*L + 2*L*K)] %*% c(alpha_0, A_mat %>% as.vector)))


# fit_FCIR_ridgelambda <- estimate_penalized_FCIR(Y = Y[,,1],
#                                                  X = X,
#                                                  Tr = Tr,
#                                                  alpha = 1e-12,
#                                                  lambda = 0.5 / (Ns * Ps),
#                                                  use_cv = FALSE)
# 
# fit_FCIR_lassolambda <- estimate_penalized_FCIR(Y = Y[,,1],
#                                                 X = X,
#                                                 Tr = Tr,
#                                                 alpha = 1,
#                                                 lambda = 0.5 / (Ns * Ps),
#                                                 use_cv = FALSE)
# 
# fit_FCIR_ridgelambdacv <- estimate_penalized_FCIR(Y = Y[,,1],
#                                                   X = X,
#                                                   Tr = Tr,
#                                                   alpha = 0,
#                                                   lambda = "lambda.min",
#                                                   cv_group_by_site = TRUE,
#                                                   use_cv = TRUE)

fit_FCIR_lassolambdacv <- estimate_penalized_FCIR(Y = Y[,,1],
                                                X = X,
                                                Tr = Tr,
                                                alpha = 1,
                                                lambda = "lambda.min",
                                                cv_group_by_site = TRUE,
                                                use_cv = TRUE)

fit_FCIR_adaptivelasso <- estimate_adaptive_lasso_FCIR(Y = Y[,,1],
                                                       X = X,
                                                       Tr,
                                                       gamma = 1, #rep(c(2,1), c(L-1 + L*K, L + L*K)),
                                                       init = "unpenalized",
                                                       lambda = "lambda.min",
                                                       use_cv = TRUE,
                                                       cv_group_by_site = TRUE,
                                                       init_lambda = "lambda.min") 

fit_FCIR_adaptivelasso2 <- estimate_adaptive_lasso_FCIR(Y = Y[,,1],
                                                       X = X,
                                                       Tr,
                                                       gamma = rep(c(2,1), c(L-1 + L*K, L + L*K)),
                                                       init = "unpenalized",
                                                       lambda = "lambda.min",
                                                       use_cv = TRUE,
                                                       cv_group_by_site = TRUE,
                                                       init_lambda = "lambda.min") 

# fit_FCIR_stepAIC <- glm(y ~ . -1, 
#                        data = data.frame(y = fit_FCIR_unpen$glm_model$y, 
#                                          fit_FCIR_unpen$glm_model %>% model.matrix),
#                        family = binomial)
# fit_FCIR_stepAIC <- MASS::stepAIC(fit_FCIR_stepAIC, direction = "backward", trace = TRUE) 
# all_names <- names(coef(fit_FCIR_unpen$glm_model))
# final_coefs <- setNames(rep(0, length(all_names)), all_names)
# final_coefs[names(coef(fit_FCIR_stepAIC))] <- coef(fit_FCIR_stepAIC)
# final_coefs


##------------------------
#' # Fit Assess results
##------------------------
beta_0; fit_FCIR_unpen$beta_0; fit_FCIR_lassolambdacv$beta_0; fit_FCIR_adaptivelasso$beta_0; fit_FCIR_adaptivelasso2$beta_0
B_mat; fit_FCIR_unpen$B_mat; fit_FCIR_lassolambdacv$B_mat; fit_FCIR_adaptivelasso$B_mat; fit_FCIR_adaptivelasso2$B_mat

alpha_0; fit_FCIR_unpen$alpha_0; fit_FCIR_lassolambdacv$alpha_0; fit_FCIR_adaptivelasso$alpha_0; fit_FCIR_adaptivelasso2$alpha_0
A_mat; fit_FCIR_unpen$A_mat; fit_FCIR_lassolambdacv$A_mat; fit_FCIR_adaptivelasso$A_mat; fit_FCIR_adaptivelasso2$A_mat




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


fit_FCIR_lassolambdacv$cv_model %>% plot
fit_FCIR_ridgelambdacv$cv_model %>% plot
fit_FCIR_lassolambdacv$fold_id

 
# gof <- -2*colSums(dbinom(fit_FCIR_unpen$glm_model$y, 
#        size = 1,
#        prob = plogis(as.matrix(model.matrix(fit_FCIR_unpen$glm_model) %*% coef(fit_FCIR_lassolambdacv$penalized_model$glmnet.fit))),
#        log = TRUE))
# 
# 
# 
# data.frame(lambda = fit_FCIR_lassolambdacv$cv_model$lambda,
#            gof = gof,
#            df = fit_FCIR_lassolambdacv$penalized_model$glmnet.fit$df,
#            AIC = gof + 2*fit_FCIR_lassolambdacv$penalized_model$glmnet.fit$df,
#            BIC = gof + log(Ns*Ps)*fit_FCIR_lassolambdacv$penalized_model$glmnet.fit$df) 
# 
