# Fit the joint Ising adaptive lasso on FCIR_M-generated Y.
# True interactions are Theta_int (site-independent); Ising Theta targets
# the same object. Species thresholds are estimated but not the comparison
# target, because FCIR_M main effects vary by site.

source("F:/ising model thesis/-Fourth-Corner-Ising-Model-/Code/Ising_joint/ising_joint_config.r")
source(file.path(ising_joint_code_dir, "estimate_Ising_joint.r"))
source(file.path(ising_joint_code_dir, "estimate_adaptive_lasso_Ising_joint.r"))

fcir_m_rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_M", "Rdata")
L = 1
est_tag_on_fcir_m = paste0("_joint_on_FCIR_M_L", L)
if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

# Match the FCIR_M simulation grid (ising_config already uses the same Ns/Ps).
Ns = c(50, 100, 200, 400, 800)
Ps = c(10, 20)

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
    data_filename = file.path(fcir_m_rdata_dir, paste0("FCIR_M_data_L", L, "_N", n, "_P", p, ".Rdata"))
    est_filename = file.path(
      rdata_dir, paste0("Ising_estimates_adaptive_lasso", est_tag_on_fcir_m, "_",
                        init_method_joint, "_N", n, "_P", p, ".Rdata")
    )

    if (!file.exists(data_filename)) {
      print(paste("FCIR_M data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste("Joint Ising-on-FCIR_M estimates already exist for N =",
                  n, ", P =", p, "- Overwriting."))
    }

    print(paste("Loading FCIR_M data for joint Ising ( N =", n, ", P =", p, ")..."))
    e = new.env()
    load(data_filename, envir = e)

    B_here = dim(e$Y)[3]
    if (nrow(e$Y) != n || ncol(e$Y) != p) {
      stop(paste0("FCIR_M data dimensions do not match N=", n, ", P=", p,
                  ": got ", paste(dim(e$Y), collapse = " x ")))
    }

    print(paste("Starting", B_here, "joint Ising adaptive lasso fits on FCIR_M data for N =",
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

    est_theta_jj = matrix(NA, nrow = B_here, ncol = p)
    est_Theta = array(NA, dim = c(p, p, B_here))
    est_lambda = numeric(B_here)

    for (b in 1:B_here) {
      est_theta_jj[b, ] = ad_results[[b]]$theta_jj
      est_Theta[, , b] = ad_results[[b]]$Theta
      est_lambda[b] = ad_results[[b]]$lambda
    }

    Theta_int = e$Theta_int
    N = n
    P = p
    B_reps = B_here
    seed = e$seed
    penalty = "adaptive_lasso"
    alpha_value = 1
    pl_form = "joint_upper_triangle"
    data_source = "FCIR_M"
    init_method = init_method_joint
    lambda_rule = lambda_rule_joint
    gamma_value = gamma_value_joint

    save(est_theta_jj, est_Theta, est_lambda,
         Theta_int,
         N, P, B_reps, seed,
         penalty, alpha_value, gamma_value, init_method, lambda_rule,
         pl_form, data_source,
         file = est_filename)
    print(paste("Saved:", est_filename))
    rm(e)
  }
}

stopCluster(cl)
print(paste("Total joint Ising-on-FCIR_M adaptive lasso execution time:",
            Sys.time() - total_start))
