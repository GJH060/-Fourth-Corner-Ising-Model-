project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_m_code_dir = file.path(project_root, "Code", "FCIR_M")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_M", "Rdata")

source(file.path(fcir_m_code_dir, "estimate_FCIR_M.r"))
source(file.path(fcir_m_code_dir, "estimate_adaptive_lasso_FCIR_M.r"))

library(foreach)
library(doParallel)

start_cluster = function(n, p) {
  n_cores = max(1, parallel::detectCores() - 1)

  print(paste("Using", n_cores, "workers for N =", n, "P =", p))
  cl = parallel::makeCluster(n_cores)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, compiler::enableJIT(0))
  parallel::clusterCall(cl, function(dir) {
    source(file.path(dir, "estimate_FCIR_M.r"))
    source(file.path(dir, "estimate_adaptive_lasso_FCIR_M.r"))
  }, dir = fcir_m_code_dir)
  cl
}

Ns = c(50, 100, 200, 400, 800)
Ps = c(10, 20)

L = 1
K = 2
B_reps = 1000
seed = 2026

gamma_value = 1
init_method = "unpenalized"
lambda_rule = "lambda.min"
unpenalize_B = TRUE
est_tag = paste0("L", L, "_unpenB")

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("FCIR_M_data_L", L, "_N", n, "_P", p, ".Rdata"))
    est_filename = file.path(rdata_dir, paste0("New_FCIR_M_estimates_adaptive_lasso_", est_tag, "_", init_method, "_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("FCIR_M data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste("FCIR_M adaptive lasso estimates already exist for N =", n, ", P =", p, "- Skipping."))
      next
    }

    print(paste("Loading FCIR_M data ( N =", n, ", P =", p, ")..."))
    load(data_filename)   # Y, X, Tr, beta_0, B_mat, Theta_int, ...

    class_counts = t(vapply(seq_len(dim(Y)[3]), function(b) {
      c(zeros = sum(Y[, , b] == 0), ones = sum(Y[, , b] == 1))
    }, numeric(2)))
    bad_reps = which(class_counts[, "zeros"] <= 1 | class_counts[, "ones"] <= 1)
    if (length(bad_reps) > 0) {
      stop(paste0(
        "FCIR_M data are too class-imbalanced for binomial glmnet at N=", n, ", P=", p,
        ". ", length(bad_reps), " replicate(s) have one class with <= 1 observation. ",
        "First bad reps: ", paste(head(bad_reps, 20), collapse = ", "),
        ". Regenerate data with less saturated Ising parameters before fitting."
      ))
    }

    print(paste("Starting", B_reps, "FCIR_M adaptive lasso fits for N =", n, "P =", p, "..."))

    est_beta_0 = matrix(NA, nrow = B_reps, ncol = L)
    est_B_mat = array(NA, dim = c(L, K, B_reps))
    est_Theta_int = array(NA, dim = c(P, P, B_reps))
    est_lambda = numeric(B_reps)

    # Ship one replicate slice per task instead of the whole Y array, so foreach
    # does not export the full N x P x B_reps array to every worker.
    Y_slices = asplit(Y, 3)
    rm(Y)

    cl = start_cluster(n, p)
    tryCatch({
      ad_results = foreach(Y_b = Y_slices, .packages = "glmnet") %dopar% {
        res = estimate_adaptive_lasso_FCIR_M(
          Y = Y_b, X = X, Tr = Tr,
          gamma = gamma_value, init = init_method,
          lambda = lambda_rule, use_cv = TRUE, cv_group_by_site = TRUE,
          unpenalize_B = unpenalize_B
        )
        list(beta_0 = res$beta_0,
             B_mat = res$B_mat,
             Theta_int = res$Theta_int,
             lambda = res$selected_lambda)
      }
    }, finally = {
      try(stopCluster(cl), silent = TRUE)
      registerDoSEQ()
    })
    rm(Y_slices)

    for (b in 1:B_reps) {
      est_beta_0[b, ] = ad_results[[b]]$beta_0
      est_B_mat[, , b] = ad_results[[b]]$B_mat
      est_Theta_int[, , b] = ad_results[[b]]$Theta_int
      est_lambda[b] = ad_results[[b]]$lambda
    }

    penalty = "adaptive_lasso"
    alpha_value = 1
    cv_type = "grouped_by_site"
    save(est_beta_0, est_B_mat, est_Theta_int, est_lambda,
         beta_0, B_mat, Theta_int,
         N = n, P = p, L, K, B_reps, seed,
         penalty, alpha_value, gamma_value, init_method, lambda_rule, cv_type,
         unpenalize_B,
         file = est_filename)
    print(paste("Saved:", est_filename))
  }
}

print(paste("Total FCIR_M adaptive lasso execution time:", Sys.time() - total_start))