project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_i_code_dir = file.path(project_root, "Code", "FCIR_I")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata_Sparse_Beta")

source(file.path(fcir_i_code_dir, "estimate_FCIR_I.r"))
source(file.path(fcir_i_code_dir, "estimate_adaptive_lasso_FCIR_I.r"))

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_FCIR_I.r"))
  source(file.path(dir, "estimate_adaptive_lasso_FCIR_I.r"))
}, dir = fcir_i_code_dir)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026

gamma_value = 1
init_method = "unpenalized"
lambda_rule = "lambda.min"

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("Rescale_FCIR_I_data_N", n, "_P", p, ".Rdata"))
    est_filename = file.path(
      rdata_dir,
      paste0("Rescale_FCIR_I_estimates_adaptive_lasso_sparse_Beta_", init_method,
             "_N", n, "_P", p, ".Rdata")
    )

    if (!file.exists(data_filename)) {
      print(paste("FCIR_I sparse-Beta data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste("FCIR_I sparse-Beta adaptive lasso estimates already exist for N =",
                  n, ", P =", p, "- Skipping."))
      next
    }

    print(paste("Loading FCIR_I sparse-Beta data ( N =", n, ", P =", p, ")..."))
    load(data_filename)   # Y, X, Tr, Beta_mat, alpha_0, A_mat, ...

    class_counts = t(vapply(seq_len(dim(Y)[3]), function(b) {
      c(zeros = sum(Y[, , b] == 0), ones = sum(Y[, , b] == 1))
    }, numeric(2)))
    bad_reps = which(class_counts[, "zeros"] <= 1 | class_counts[, "ones"] <= 1)
    if (length(bad_reps) > 0) {
      stop(paste0(
        "FCIR_I sparse-Beta data are too class-imbalanced for binomial glmnet at N=", n,
        ", P=", p, ". ", length(bad_reps),
        " replicate(s) have one class with <= 1 observation. ",
        "First bad reps: ", paste(head(bad_reps, 20), collapse = ", "),
        ". Regenerate data with less saturated Ising parameters before fitting."
      ))
    }

    print(paste("Starting", B_reps, "FCIR_I sparse-Beta adaptive lasso fits for N =",
                n, "P =", p, "..."))

    est_Beta_mat = array(NA, dim = c(P, L, B_reps))
    est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
    est_A_mat = array(NA, dim = c(L, K, B_reps))
    est_lambda = numeric(B_reps)

    Y_slices = asplit(Y, 3)
    rm(Y)

    ad_results = foreach(Y_b = Y_slices, .packages = "glmnet") %dopar% {
      res = estimate_adaptive_lasso_FCIR_I(
        Y = Y_b, X = X, Tr = Tr,
        gamma = gamma_value, init = init_method,
        lambda = lambda_rule, use_cv = TRUE, cv_group_by_site = TRUE
      )
      list(Beta_mat = res$Beta_mat,
           alpha_0 = res$alpha_0,
           A_mat = res$A_mat,
           lambda = res$selected_lambda)
    }
    rm(Y_slices)

    for (b in 1:B_reps) {
      est_Beta_mat[, , b] = ad_results[[b]]$Beta_mat
      est_alpha_0[b, ] = ad_results[[b]]$alpha_0
      est_A_mat[, , b] = ad_results[[b]]$A_mat
      est_lambda[b] = ad_results[[b]]$lambda
    }

    penalty = "adaptive_lasso"
    alpha_value = 1
    cv_type = "grouped_by_site"
    beta_structure = "sparse"
    save(est_Beta_mat, est_alpha_0, est_A_mat, est_lambda,
         Beta_mat, alpha_0, A_mat,
         N = n, P = p, L, K, B_reps, seed,
         penalty, alpha_value, gamma_value, init_method, lambda_rule, cv_type,
         beta_structure, file = est_filename)
    print(paste("Saved:", est_filename))
  }
}

stopCluster(cl)
print(paste("Total FCIR_I sparse-Beta adaptive lasso execution time:",
            Sys.time() - total_start))
