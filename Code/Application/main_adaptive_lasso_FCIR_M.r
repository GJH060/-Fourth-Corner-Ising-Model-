# Adaptive lasso MPLE of the FCIR_M model on the Southern Ocean CPR data.
# Self-contained: sets its own paths, builds its own inputs, fits, saves.
#
#   Rscript Code/Application/main_adaptive_lasso_FCIR_M.r
#   source("Code/Application/main_adaptive_lasso_FCIR_M.r")
#
# FCIR_M keeps the fourth-corner main effect but replaces the trait-driven
# interaction by a free, site-independent residual association matrix -- here
# P(P-1)/2 = 435 edges, which is where the adaptive lasso does most of its work.

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_m_code_dir = file.path(project_root, "Code", "FCIR_M")
app_data_dir = file.path(project_root, "Application")
rdata_dir = file.path(project_root, "Simulation_Results", "Application", "Rdata")
if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

source(file.path(fcir_m_code_dir, "estimate_FCIR_M.r"))
source(file.path(fcir_m_code_dir, "estimate_adaptive_lasso_FCIR_M.r"))
library(glmnet)
library(Matrix)

# --- Settings ---------------------------------------------------------------
env_vars = c("salinity", "water_temperature", "photosynthetically_active_radiation")

N_use = NA          # NA = all 29,983 sites; set an integer for a subsample
seed = 2026

gamma_value = 1
init_method = "unpenalized"
lambda_rule = "lambda.min"
cv_type = "grouped_by_site"

# The simulations used unpenalize_B = TRUE, but Traitmat here is a full
# indicator basis, so the main-effect block has an L-dimensional null space.
# Leaving it unpenalized would leave that null space completely unregularized
# and coordinate descent free to drift along it. Penalize it.
unpenalize_B = FALSE

# Build the 435 edge columns as a dgCMatrix. Dense, they are two N*P x 435
# outer() temporaries -- about 15 GB at the full sample. Sparse they hold
# sum_s R_s * (P-1) non-zeros, roughly 135 MB, and the fit is numerically
# identical (verified: design matrices equal elementwise, lambda.min equal to
# 8 digits, coefficients agreeing to 5e-14).
use_sparse = TRUE

# --- Data -------------------------------------------------------------------
load(file.path(app_data_dir, "reduceddat.RData"))   # resp, X, Traitmat

# Y must be a numeric matrix, and the first column of X must be constant --
# see the long note in main_adaptive_lasso_FCIR.r.
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
n_edges = P * (P - 1) / 2
n_params = L + L * K + n_edges

est_filename = file.path(rdata_dir,
                         paste0("app_FCIR_M_adaptive_lasso_", init_method, "_N", N, ".Rdata"))

if (file.exists(est_filename)) {
  print(paste("FCIR_M estimates already exist - Skipping:", est_filename))
} else {

  print(paste0("FCIR_M: N = ", N, ", P = ", P, ", L = ", L, ", K = ", K,
               " (", format(N * P, big.mark = ","), " pseudo-likelihood rows, ",
               n_params, " design columns, ", n_edges, " free edges)"))

  fit_warnings = character(0)
  total_start = Sys.time()
  ad_fit = withCallingHandlers(
    estimate_adaptive_lasso_FCIR_M(Y = Y, X = Xd, Tr = Tr,
                                   gamma = gamma_value,
                                   init = init_method,
                                   lambda = lambda_rule,
                                   use_cv = TRUE,
                                   cv_group_by_site = TRUE,
                                   unpenalize_B = unpenalize_B,
                                   sparse = use_sparse),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  elapsed = as.numeric(difftime(Sys.time(), total_start, units = "secs"))

  est_beta_0 = ad_fit$beta_0
  est_B_mat = ad_fit$B_mat
  est_Theta_int = ad_fit$Theta_int
  est_lambda = ad_fit$selected_lambda
  init_rank = ad_fit$init_model$glm_model$rank
  init_aliased = sum(is.na(coef(ad_fit$init_model$glm_model)))
  penalty_factor = ad_fit$penalty_factor

  edges = est_Theta_int[upper.tri(est_Theta_int)]
  print(paste("lambda.min =", signif(est_lambda, 5)))
  print(paste0("glm init: rank ", init_rank, " of ", n_params,
               " columns, ", init_aliased, " aliased"))
  print(paste0("edges selected: ", sum(edges != 0), "/", n_edges,
               "  (+", sum(edges > 0), " / -", sum(edges < 0), ")"))
  print(paste0("zeros: B ", sum(est_B_mat == 0), "/", L * K))
  if (length(fit_warnings) > 0) for (w in unique(fit_warnings)) print(paste("WARNING:", w))
  print(paste("Elapsed:", sprintf("%.1f", elapsed), "s"))

  penalty = "adaptive_lasso"
  alpha_value = 1
  save(est_beta_0, est_B_mat, est_Theta_int, est_lambda,
       init_rank, init_aliased, penalty_factor, fit_warnings,
       N, P, L, K, n_edges, env_vars, N_use, seed, elapsed,
       penalty, alpha_value, gamma_value, init_method, lambda_rule, cv_type,
       unpenalize_B, use_sparse,
       file = est_filename)
  print(paste("Saved:", est_filename))
}
