source("F:/ising model thesis/-Fourth-Corner-Ising-Model-/Code/Ising/ising_config.r")
source(file.path(ising_code_dir, "estimate_Ising.r"))
source(file.path(ising_code_dir, "estimate_adaptive_lasso_Ising.r"))

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_Ising.r"))
  source(file.path(dir, "estimate_adaptive_lasso_Ising.r"))
}, dir = ising_code_dir)

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(
      rdata_dir, paste0("Ising_data", setting_tag, "_N", n, "_P", p, ".Rdata")
    )
    est_filename = file.path(
      rdata_dir, paste0("Ising_estimates_adaptive_lasso", setting_tag, "_",
                        init_method, "_N", n, "_P", p, ".Rdata")
    )

    if (!file.exists(data_filename)) {
      print(paste("Ising data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste("Ising adaptive lasso estimates already exist for N =", n,
                  ", P =", p, "- Skipping."))
      next
    }

    print(paste("Loading Ising data ( N =", n, ", P =", p, ")..."))
    e = new.env()
    load(data_filename, envir = e)   # Y, theta_jj, Theta, ...

    print(paste("Starting", B_reps, "Ising adaptive lasso fits for N =", n,
                "P =", p, "..."))

    Y_slices = asplit(e$Y, 3)

    ad_results = foreach(Y_b = Y_slices, .packages = "glmnet") %dopar% {
      res = estimate_adaptive_lasso_Ising(
        Y = Y_b,
        gamma = gamma_value, init = init_method,
        lambda = lambda_rule, use_cv = TRUE,
        symmetrize = symmetrize_rule
      )
      list(theta_jj = res$theta_jj,
           Theta = res$Theta,
           lambda = res$selected_lambda,
           n_degenerate = sum(res$degenerate_nodes))
    }
    rm(Y_slices)

    est_theta_jj = matrix(NA, nrow = B_reps, ncol = p)
    est_Theta = array(NA, dim = c(p, p, B_reps))
    est_lambda = matrix(NA, nrow = B_reps, ncol = p)
    est_n_degenerate = numeric(B_reps)

    for (b in 1:B_reps) {
      est_theta_jj[b, ] = ad_results[[b]]$theta_jj
      est_Theta[, , b] = ad_results[[b]]$Theta
      est_lambda[b, ] = ad_results[[b]]$lambda
      est_n_degenerate[b] = ad_results[[b]]$n_degenerate
    }

    theta_jj = e$theta_jj
    Theta = e$Theta
    N = n
    P = p
    penalty = "adaptive_lasso"
    alpha_value = 1

    save(est_theta_jj, est_Theta, est_lambda, est_n_degenerate,
         theta_jj, Theta,
         N, P, B_reps, seed,
         edge_density, edge_values, edge_magnitude,
         penalty, alpha_value, gamma_value, init_method, lambda_rule,
         symmetrize_rule,
         file = est_filename)
    n_bad = sum(est_n_degenerate)
    if (n_bad > 0) {
      print(paste0("  NOTE: ", n_bad, " of ", B_reps * p,
                   " node fits were degenerate (too few observations in one",
                   " class to cross-validate) and contributed no edges; ",
                   sum(est_n_degenerate > 0), " of ", B_reps,
                   " replicates affected."))
    }
    print(paste("Saved:", est_filename))
    rm(e)
  }
}

stopCluster(cl)
print(paste("Total Ising adaptive lasso execution time:",
            Sys.time() - total_start))
