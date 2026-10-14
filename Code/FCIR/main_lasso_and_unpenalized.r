# Comparison estimators for the FCIR simulation, fitted to the same data as
# main_adaptive_lasso.r:
#   - ordinary lasso (unit penalty weights, otherwise identical pipeline)
#   - unpenalized MPLE (plain stacked logistic regression, original scale)

project_root = "D:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_code_dir = file.path(project_root, "Code", "FCIR")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Sparse")

source(file.path(fcir_code_dir, "estimate_FCIR.r"))
source(file.path(fcir_code_dir, "estimate_adaptive_lasso_FCIR.r"))

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
# The byte-code JIT sporadically fails inside the design-building loops with
# "cannot change value of locked binding"; disable it on the workers.
clusterEvalQ(cl, compiler::enableJIT(0))
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_FCIR.r"))
  source(file.path(dir, "estimate_adaptive_lasso_FCIR.r"))
}, dir = fcir_code_dir)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026
lambda_rule = "lambda.min"

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("FCIR_data_N", n, "_P", p, ".Rdata"))
    lasso_filename = file.path(rdata_dir, paste0("FCIR_estimates_lasso_N", n, "_P", p, ".Rdata"))
    unpen_filename = file.path(rdata_dir, paste0("FCIR_estimates_unpenalized_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("Data not found, skipping:", data_filename))
      next
    }
    if (file.exists(lasso_filename) && file.exists(unpen_filename)) {
      print(paste("Lasso and unpenalized estimates already exist for N =", n, ", P =", p, "- Skipping."))
      next
    }

    load(data_filename)   # Y, X, Tr, beta_0, B_mat, alpha_0, A_mat, ...
    print(paste("Starting", B_reps, "lasso + unpenalized fits for N =", n, "P =", p, "..."))

    Y_slices = asplit(Y, 3)
    fits = foreach(b = seq_len(B_reps), Y_b = Y_slices, .packages = "glmnet") %dopar% {
      # Retry a replicate on a sporadic native/JIT error instead of losing the
      # whole configuration; reseeding keeps the CV folds identical.
      for (attempt in 1:3) {
        out = tryCatch({
          set.seed(seed + b)   # reproducible CV folds
          list(unpen = suppressWarnings(estimate_unpenalized_FCIR(Y = Y_b, X = X, Tr = Tr, standardize = FALSE)),
               las = suppressWarnings(estimate_adaptive_lasso_FCIR(
                 Y = Y_b, X = X, Tr = Tr,
                 lambda = lambda_rule, use_cv = TRUE, cv_group_by_site = TRUE,
                 adaptive = FALSE)))
        }, error = function(e) if (attempt == 3) stop(e) else NULL)
        if (!is.null(out)) break
      }
      unpen = out$unpen
      las = out$las
      list(unpen = list(beta_0 = unname(unpen$beta_0), B_mat = unpen$B_mat,
                        alpha_0 = unname(unpen$alpha_0), A_mat = unpen$A_mat,
                        converged = unpen$glm_model$converged),
           lasso = list(beta_0 = las$beta_0, B_mat = las$B_mat,
                        alpha_0 = las$alpha_0, A_mat = las$A_mat,
                        lambda = las$selected_lambda))
    }
    rm(Y_slices)

    collect = function(key) {
      list(
        est_beta_0  = t(vapply(fits, function(f) f[[key]]$beta_0, numeric(L))),
        est_B_mat   = simplify2array(lapply(fits, function(f) f[[key]]$B_mat)),
        est_alpha_0 = t(vapply(fits, function(f) f[[key]]$alpha_0, numeric(L))),
        est_A_mat   = simplify2array(lapply(fits, function(f) f[[key]]$A_mat))
      )
    }

    with(collect("lasso"), {
      est_lambda = vapply(fits, function(f) f$lasso$lambda, numeric(1))
      penalty = "lasso"; alpha_value = 1; cv_type = "grouped_by_site"
      save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat, est_lambda,
           beta_0, B_mat, alpha_0, A_mat,
           N = n, P = p, L, K, B_reps, seed, penalty, alpha_value, lambda_rule, cv_type,
           file = lasso_filename)
    })
    with(collect("unpen"), {
      est_converged = vapply(fits, function(f) f$unpen$converged, logical(1))
      method = "unpenalized"
      save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat, est_converged,
           beta_0, B_mat, alpha_0, A_mat,
           N = n, P = p, L, K, B_reps, seed, method,
           file = unpen_filename)
    })
    rm(fits)
    print(paste("Saved:", lasso_filename, "and", unpen_filename))
  }
}

stopCluster(cl)
print(paste("Total lasso + unpenalized execution time:", format(Sys.time() - total_start)))
