# Adaptive lasso MPLE of the FCIR model on the Southern Ocean CPR data.
# Self-contained: sets its own paths, builds its own inputs, fits, saves.
#
#   Rscript Code/Application/main_adaptive_lasso_FCIR.r
#   source("Code/Application/main_adaptive_lasso_FCIR.r")

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_code_dir = file.path(project_root, "Code", "FCIR")
app_data_dir = file.path(project_root, "Application")
rdata_dir = file.path(project_root, "Simulation_Results", "Application", "Rdata")
if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

source(file.path(fcir_code_dir, "estimate_FCIR.r"))
source(file.path(fcir_code_dir, "estimate_adaptive_lasso_FCIR.r"))
library(glmnet)

# --- Settings ---------------------------------------------------------------
# x_s = (1, salinity, water_temperature, PAR)  ->  L = 4.
# Latitude/longitude are excluded: space-time correlation is out of scope here,
# and water_temperature already correlates 0.87 with latitude.
env_vars = c("salinity", "water_temperature", "photosynthetically_active_radiation")

N_use = NA          # NA = all 29,983 sites; set an integer for a subsample
seed = 2026

gamma_value = 1
init_method = "unpenalized"
lambda_rule = "lambda.min"
cv_type = "grouped_by_site"

# --- Data -------------------------------------------------------------------
load(file.path(app_data_dir, "reduceddat.RData"))   # resp, X, Traitmat

# Two things estimate_FCIR.r assumes but does not check.
#
# (1) Y must be a numeric matrix. It indexes rows as y_s = Y[s, ] and feeds x_s
#     to kronecker(); a data.frame row is itself a data.frame and kronecker()
#     fails with "invalid 'dimnames' given for data frame".
#
# (2) The first column of X must be constant. estimate_FCIR.r fits
#         glm(y ~ ., data = data.frame(y = glm_Y, glm_X[, -1]))
#     which unconditionally drops design column 1 and lets glm supply its own
#     intercept -- only equivalent to FCIR when column 1 IS the ones column.
#     Pass three bare covariates and it returns the intercept in beta_0[1],
#     never estimates the first covariate, and shifts every parameter slice by
#     one, with no error and no warning. The stopifnot below is what catches it.
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

est_filename = file.path(rdata_dir,
                         paste0("app_FCIR_adaptive_lasso_", init_method, "_N", N, ".Rdata"))

if (file.exists(est_filename)) {
  print(paste("FCIR estimates already exist - Skipping:", est_filename))
} else {

  print(paste0("FCIR: N = ", N, ", P = ", P, ", L = ", L, ", K = ", K,
               " (", format(N * P, big.mark = ","), " pseudo-likelihood rows, ",
               2 * L + 2 * L * K, " design columns)"))

  # Traitmat is the full indicator expansion of broad_func_grp, so the K trait
  # blocks of the design sum to the beta_0 block and the main-effect block is
  # rank deficient by L. glm()'s QR aliases the last trait block, whose adaptive
  # weight then becomes 1/eps and forces B[, K] to exactly zero -- i.e. the last
  # trait column acts as the reference group. Expect 4 aliased columns.
  # Quasi-separation ("fitted probabilities numerically 0 or 1") is the warning
  # that actually matters: it drives adaptive weights toward zero and lets those
  # columns escape the penalty.
  fit_warnings = character(0)
  total_start = Sys.time()
  ad_fit = withCallingHandlers(
    estimate_adaptive_lasso_FCIR(Y = Y, X = Xd, Tr = Tr,
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

  est_beta_0 = ad_fit$beta_0
  est_B_mat = ad_fit$B_mat
  est_alpha_0 = ad_fit$alpha_0
  est_A_mat = ad_fit$A_mat
  est_lambda = ad_fit$selected_lambda
  init_rank = ad_fit$init_model$glm_model$rank
  init_aliased = sum(is.na(coef(ad_fit$init_model$glm_model)))
  penalty_factor = ad_fit$penalty_factor

  print(paste("lambda.min =", signif(est_lambda, 5)))
  print(paste0("glm init: rank ", init_rank, " of ", 2 * L + 2 * L * K,
               " columns, ", init_aliased, " aliased"))
  print(paste0("zeros: B ", sum(est_B_mat == 0), "/", L * K,
               ", A ", sum(est_A_mat == 0), "/", L * K,
               ", alpha_0 ", sum(est_alpha_0 == 0), "/", L))
  if (length(fit_warnings) > 0) for (w in unique(fit_warnings)) print(paste("WARNING:", w))
  print(paste("Elapsed:", sprintf("%.1f", elapsed), "s"))

  # Save the parameter blocks only. ad_fit$init_model is a full glm object
  # carrying an NP x n_par QR factor -- hundreds of MB to several GB -- while
  # the estimates themselves are a couple of KB.
  penalty = "adaptive_lasso"
  alpha_value = 1
  save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat, est_lambda,
       init_rank, init_aliased, penalty_factor, fit_warnings,
       N, P, L, K, env_vars, N_use, seed, elapsed,
       penalty, alpha_value, gamma_value, init_method, lambda_rule, cv_type,
       file = est_filename)
  print(paste("Saved:", est_filename))
}
