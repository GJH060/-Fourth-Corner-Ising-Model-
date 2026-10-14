# Unpenalized MPLE for the L = 3 FCIR_M simulation (data from
# main_generate_data_FCIR_M_L3.r), used as the RMSE baseline for the adaptive lasso.

project_root = "D:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_m_code_dir = file.path(project_root, "Code", "FCIR_M")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_M", "Rdata")

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, compiler::enableJIT(0))
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_FCIR_M.r"))
}, dir = fcir_m_code_dir)

Ns = c(50, 100, 200, 400, 800)
Ps = c(10, 20)

L = 3
K = 2
B_reps = 1000
seed = 2026

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("FCIR_M_data_L", L, "_N", n, "_P", p, ".Rdata"))
    est_filename = file.path(rdata_dir, paste0("FCIR_M_estimates_unpenalized_L", L, "_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("FCIR_M data not found, skipping:", data_filename))
      next
    }
    if (file.exists(est_filename)) {
      print(paste("FCIR_M unpenalized estimates already exist for N =", n, ", P =", p, "- Skipping."))
      next
    }

    load(data_filename)   # Y, X, Tr, beta_0, B_mat, Theta_int, ...
    print(paste("Starting", B_reps, "FCIR_M unpenalized fits for N =", n, "P =", p, "..."))

    Y_slices = asplit(Y, 3)
    fits = foreach(Y_b = Y_slices, .export = c("X", "Tr")) %dopar% {
      res = suppressWarnings(estimate_unpenalized_FCIR_M(Y = Y_b, X = X, Tr = Tr))
      coefs = res$glm_model$coefficients
      list(beta_0 = unname(res$beta_0), B_mat = res$B_mat, Theta_int = res$Theta_int,
           converged = res$glm_model$converged, n_aliased = sum(is.na(coefs)))
    }
    rm(Y_slices)

    est_beta_0 = t(vapply(fits, function(f) f$beta_0, numeric(L)))
    est_B_mat = simplify2array(lapply(fits, function(f) f$B_mat))
    est_Theta_int = simplify2array(lapply(fits, function(f) f$Theta_int))
    est_converged = vapply(fits, function(f) f$converged, logical(1))
    est_n_aliased = vapply(fits, function(f) f$n_aliased, numeric(1))
    rm(fits)

    method = "unpenalized"
    save(est_beta_0, est_B_mat, est_Theta_int, est_converged, est_n_aliased,
         beta_0, B_mat, Theta_int,
         N = n, P = p, L, K, B_reps, seed, method,
         file = est_filename)
    print(paste("Saved:", est_filename, "- converged:", sum(est_converged), "/", B_reps,
                "- reps with aliased edges:", sum(est_n_aliased > 0)))
  }
}

stopCluster(cl)
print(paste("Total FCIR_M unpenalized execution time:", format(Sys.time() - total_start)))
