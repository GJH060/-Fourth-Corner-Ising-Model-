project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_code_dir = file.path(project_root, "Code", "FCIR")

use_dense = FALSE   # FALSE -> Rdata_Sparse, TRUE -> Rdata_Dense

rdata_dir = file.path(project_root, "Simulation_Results",
                      if (use_dense) "Rdata_Dense" else "Rdata_Sparse")
infix = if (use_dense) "_Dense" else ""

source(file.path(fcir_code_dir, "estimate_adaptive_lasso_FCIR.r"))

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_adaptive_lasso_FCIR.r"))
}, dir = fcir_code_dir)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026

# Adaptive lasso settings
gamma_value = 1
init_method = "unpenalized"
lambda_rule = "lambda.min"

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {

    data_filename = file.path(rdata_dir, paste0("FCIR_data", infix, "_N", n, "_P", p, ".Rdata"))
    est_filename  = file.path(rdata_dir, paste0("FCIR_estimates_adaptive_lasso_", init_method, infix, "_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("Data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste("Adaptive lasso estimates already exist for N =", n, ", P =", p, "- Skipping."))
      next
    }

    print(paste("Loading data ( N =", n, ", P =", p, ")..."))
    load(data_filename)   # Y, X, Tr, beta_0, B_mat, alpha_0, A_mat, ...

    print(paste("Starting", B_reps, "adaptive lasso fits for N =", n, "P =", p, "..."))

    est_beta_0  = matrix(NA, nrow = B_reps, ncol = L)
    est_B_mat   = array(NA, dim = c(L, K, B_reps))
    est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
    est_A_mat   = array(NA, dim = c(L, K, B_reps))
    est_lambda  = numeric(B_reps)

    ad_results = foreach(b = 1:B_reps, .packages = "glmnet") %dopar% {
      res = estimate_adaptive_lasso_FCIR(
        Y = Y[, , b], X = X, Tr = Tr,
        gamma = gamma_value, init = init_method,
        lambda = lambda_rule, use_cv = TRUE, cv_group_by_site = TRUE
      )
      list(beta_0 = res$beta_0, B_mat = res$B_mat,
           alpha_0 = res$alpha_0, A_mat = res$A_mat,
           lambda = res$selected_lambda)
    }

    for (b in 1:B_reps) {
      est_beta_0[b, ]  = ad_results[[b]]$beta_0
      est_B_mat[, , b] = ad_results[[b]]$B_mat
      est_alpha_0[b, ] = ad_results[[b]]$alpha_0
      est_A_mat[, , b] = ad_results[[b]]$A_mat
      est_lambda[b]    = ad_results[[b]]$lambda
    }

    penalty = "adaptive_lasso"
    alpha_value = 1
    cv_type = "grouped_by_site"
    save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat, est_lambda,
         beta_0, B_mat, alpha_0, A_mat,
         N = n, P = p, L, K, B_reps, seed,
         penalty, alpha_value, gamma_value, init_method, lambda_rule, cv_type,
         file = est_filename)
    print(paste("Saved:", est_filename))
  }
}

stopCluster(cl)
print(paste("Total adaptive lasso execution time:", Sys.time() - total_start))