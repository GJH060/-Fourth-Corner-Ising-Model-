source("F:/ising model thesis/-Fourth-Corner-Ising-Model-/Code/Ising_joint/ising_joint_config.r")
source(file.path(ising_joint_code_dir, "estimate_Ising_joint.r"))
source(file.path(ising_joint_code_dir, "estimate_adaptive_lasso_Ising_joint.r"))

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_Ising_joint.r"))
  source(file.path(dir, "estimate_adaptive_lasso_Ising_joint.r"))
}, dir = ising_joint_code_dir)

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    # Same generated data as the node-wise Ising pipeline.
    data_filename = file.path(
      rdata_dir, paste0("Ising_data", setting_tag, "_N", n, "_P", p, ".Rdata")
    )
    est_filename = file.path(
      rdata_dir, paste0("Ising_estimates_adaptive_lasso", est_tag_joint, "_",
                        init_method_joint, "_N", n, "_P", p, ".Rdata")
    )

    if (!file.exists(data_filename)) {
      print(paste("Ising data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste("Joint Ising adaptive lasso estimates already exist for N =",
                  n, ", P =", p, "- Skipping."))
      next
    }

    print(paste("Loading Ising data ( N =", n, ", P =", p, ")..."))
    e = new.env()
    load(data_filename, envir = e)

    print(paste("Starting", B_reps, "joint Ising adaptive lasso fits for N =",
                n, "P =", p, "..."))

    Y_slices = asplit(e$Y, 3)

    ad_results = foreach(Y_b = Y_slices, .packages = "glmnet") %dopar% {
      res = estimate_adaptive_lasso_Ising_joint(
        Y = Y_b,
        gamma = gamma_value_joint,
        init = init_method_joint,
        lambda = lambda_rule_joint,
        use_cv = TRUE,
        cv_group_by_site = TRUE
      )
      list(theta_jj = res$theta_jj,
           Theta = res$Theta,
           lambda = res$selected_lambda)
    }
    rm(Y_slices)

    est_theta_jj = matrix(NA, nrow = B_reps, ncol = p)
    est_Theta = array(NA, dim = c(p, p, B_reps))
    est_lambda = numeric(B_reps)

    for (b in 1:B_reps) {
      est_theta_jj[b, ] = ad_results[[b]]$theta_jj
      est_Theta[, , b] = ad_results[[b]]$Theta
      est_lambda[b] = ad_results[[b]]$lambda
    }

    theta_jj = e$theta_jj
    Theta = e$Theta
    N = n
    P = p
    penalty = "adaptive_lasso"
    alpha_value = 1
    pl_form = "joint_upper_triangle"
    init_method = init_method_joint
    lambda_rule = lambda_rule_joint
    gamma_value = gamma_value_joint

    save(est_theta_jj, est_Theta, est_lambda,
         theta_jj, Theta,
         N, P, B_reps, seed,
         edge_density, edge_values, edge_magnitude,
         penalty, alpha_value, gamma_value, init_method, lambda_rule,
         pl_form,
         file = est_filename)
    print(paste("Saved:", est_filename))
    rm(e)
  }
}

stopCluster(cl)
print(paste("Total joint Ising adaptive lasso execution time:",
            Sys.time() - total_start))
