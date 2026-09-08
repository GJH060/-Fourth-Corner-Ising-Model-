# Adaptive lasso MPLE of the FCIR_I model on the Southern Ocean CPR data.
# Self-contained: sets its own paths, builds its own inputs, fits, saves.
#
#   Rscript Code/Application/main_adaptive_lasso_FCIR_I.r
#   source("Code/Application/main_adaptive_lasso_FCIR_I.r")
#
# FCIR_I replaces the fourth-corner main effect by unconstrained species-specific
# environmental responses beta_j, keeping the trait-driven interaction. It has no
# beta_0, so unlike FCIR and FCIR_M its design is full rank and no trait column
# acts as a reference group.

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_i_code_dir = file.path(project_root, "Code", "FCIR_I")
app_data_dir = file.path(project_root, "Application")
rdata_dir = file.path(project_root, "Simulation_Results", "Application", "Rdata")
if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

source(file.path(fcir_i_code_dir, "estimate_FCIR_I.r"))
source(file.path(fcir_i_code_dir, "estimate_adaptive_lasso_FCIR_I.r"))
library(glmnet)

# --- Settings ---------------------------------------------------------------
env_vars = c("salinity", "water_temperature", "photosynthetically_active_radiation")

N_use = NA          # NA = all 29,983 sites; set an integer for a subsample
seed = 2026

gamma_value = 1
init_method = "unpenalized"
lambda_rule = "lambda.min"
cv_type = "grouped_by_site"

# --- Data -------------------------------------------------------------------
load(file.path(app_data_dir, "reduceddat.RData"))   # resp, X, Traitmat

# Y must be a numeric matrix, and the first column of X must be constant --
# see the long note in main_adaptive_lasso_FCIR.r. estimate_FCIR_I.r fits
# glm(glm_Y ~ glm_X + 0) with no separate intercept, so requirement (2) is not
# the same silent trap it is for FCIR, but x_s still has to carry the ones
# column for beta_j to contain a species intercept at all.
Y = as.matrix(resp)
storage.mode(Y) = "numeric"
Xd = cbind(1, as.matrix(X[, env_vars]))
colnames(Xd) = c("(const)", env_vars)
Tr = Traitmat

stopifnot(is.matrix(Y), is.numeric(Y), all(Y %in% c(0, 1)),
          is.matrix(Xd), is.numeric(Xd),
          sd(Xd[, 1]) == 0, Xd[1, 1] == 1,
          nrow(Y) == nrow(Xd),
          identical(colnames(Y), rownames(Tr)))

if (!is.na(N_use)) {
  set.seed(seed)
  ix = sort(sample(nrow(Y), N_use))
  Y = Y[ix, , drop = FALSE]
  Xd = Xd[ix, , drop = FALSE]
}

N = nrow(Y); P = ncol(Y); L = ncol(Xd); K = ncol(Tr)
n_params = P * L + L + L * K

est_filename = file.path(rdata_dir,
                         paste0("app_FCIR_I_adaptive_lasso_", init_method, "_N", N, ".Rdata"))

if (file.exists(est_filename)) {
  print(paste("FCIR_I estimates already exist - Skipping:", est_filename))
} else {

  print(paste0("FCIR_I: N = ", N, ", P = ", P, ", L = ", L, ", K = ", K,
               " (", format(N * P, big.mark = ","), " pseudo-likelihood rows, ",
               n_params, " design columns)"))

  fit_warnings = character(0)
  total_start = Sys.time()
  ad_fit = withCallingHandlers(
    estimate_adaptive_lasso_FCIR_I(Y = Y, X = Xd, Tr = Tr,
                                   gamma = gamma_value,
                                   init = init_method,
                                   lambda = lambda_rule,
                                   use_cv = TRUE,
                                   cv_group_by_site = TRUE),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  elapsed = as.numeric(difftime(Sys.time(), total_start, units = "secs"))

  est_Beta_mat = ad_fit$Beta_mat          # P x L, row j is beta_j
  est_alpha_0 = ad_fit$alpha_0
  est_A_mat = ad_fit$A_mat
  est_lambda = ad_fit$selected_lambda
  init_rank = ad_fit$init_model$glm_model$rank
  init_aliased = sum(is.na(coef(ad_fit$init_model$glm_model)))
  penalty_factor = ad_fit$penalty_factor

  print(paste("lambda.min =", signif(est_lambda, 5)))
  print(paste0("glm init: rank ", init_rank, " of ", n_params,
               " columns, ", init_aliased, " aliased (expect 0 for FCIR_I)"))
  print(paste0("zeros: Beta ", sum(est_Beta_mat == 0), "/", P * L,
               ", A ", sum(est_A_mat == 0), "/", L * K,
               ", alpha_0 ", sum(est_alpha_0 == 0), "/", L))
  if (length(fit_warnings) > 0) for (w in unique(fit_warnings)) print(paste("WARNING:", w))
  print(paste("Elapsed:", sprintf("%.1f", elapsed), "s"))

  penalty = "adaptive_lasso"
  alpha_value = 1
  save(est_Beta_mat, est_alpha_0, est_A_mat, est_lambda,
       init_rank, init_aliased, penalty_factor, fit_warnings,
       N, P, L, K, env_vars, N_use, seed, elapsed,
       penalty, alpha_value, gamma_value, init_method, lambda_rule, cv_type,
       file = est_filename)
  print(paste("Saved:", est_filename))
}
