rm(list = ls())
library(tidyverse)
library(glmnet)
source("Code/FCIR/estimate_FCIR.r")

load("../application_SOCPR/reduceddat.RData")

Traitmat <- Traitmat %>%
    as.data.frame() 
functional_group <- factor(names(Traitmat)[max.col(as.matrix(Traitmat))], 
                           levels = names(Traitmat)) 
Tr <- model.matrix(~ 0 + functional.group, data = data.frame(functional.group = functional_group))[,-1]

table(X$year)
sel_year <- which(X$year %in% c(2014:2016))
Y <- as.matrix(resp[sel_year,])
Xd <- model.matrix(~ salinity + water_temperature + photosynthetically_active_radiation, data = X[sel_year,]) %>% 
    as.matrix


##------------------------
#' # Fit various forms of the FCIR model or a special case of it
##------------------------
fit_FCIR_unpen_nostan <- estimate_unpenalized_FCIR(Y = Y,
                                                   X = Xd,
                                                   Tr = Tr,
                                                   standardize = FALSE)


fit_FCIR_unpen_nostan$beta_0
fit_FCIR_unpen_nostan$B_mat
fit_FCIR_unpen_nostan$A_mat
fit_FCIR_unpen_nostan$alpha_0
