# Ordinary lasso fits for the FCIR_I simulation reported in the thesis (sparse
# Beta, rescaled data in Rdata_Sparse_Beta), as a comparison estimator for the
# adaptive lasso of main_adaptive_lasso_sparse_Beta_FCIR_I.r. Same data, same
# standardization, same site-grouped CV; the only difference is that every
# penalized coefficient gets penalty weight 1 (adaptive = FALSE). The
# species-specific intercepts stay unpenalized.
#
# Unlike the adaptive lasso script, each replicate sets seed + b before fitting,
# so the CV folds (and hence the results) are reproducible.
#
# Computation time is recorded: est_time_sec holds the elapsed seconds of each
# replicate's fit (measured inside the worker), wall_time_sec the elapsed time
# of the whole configuration with n_workers parallel workers.
#
# Set the environment variable FCIR_PROJECT_ROOT to the repository root if it
# differs from the default below.
project_root = Sys.getenv("FCIR_PROJECT_ROOT",
                          unset = "D:/ising model thesis/-Fourth-Corner-Ising-Model-")
fcir_i_code_dir = file.path(project_root, "Code", "FCIR_I")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata_Sparse_Beta")

source(file.path(fcir_i_code_dir, "estimate_FCIR_I.r"))
source(file.path(fcir_i_code_dir, "estimate_adaptive_lasso_FCIR_I.r"))

library(foreach)
library(doParallel)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026
lambda_rule = "lambda.min"

# Optional: restrict the run, e.g. FCIR_I_NS="800" FCIR_I_PS="30" Rscript main_lasso_sparse_Beta_FCIR_I.r
if (nzchar(Sys.getenv("FCIR_I_NS"))) Ns = as.numeric(strsplit(Sys.getenv("FCIR_I_NS"), ",")[[1]])
if (nzchar(Sys.getenv("FCIR_I_PS"))) Ps = as.numeric(strsplit(Sys.getenv("FCIR_I_PS"), ",")[[1]])
if (nzchar(Sys.getenv("FCIR_B_REPS"))) B_reps = as.integer(Sys.getenv("FCIR_B_REPS"))

n_workers = max(1, parallel::detectCores() - 1)
machine_info = list(cpu = Sys.getenv("PROCESSOR_IDENTIFIER"),
                    logical_cores = parallel::detectCores(),
                    os = paste(Sys.info()[c("sysname", "release")], collapse = " "),
                    r_version = R.version.string,
                    glmnet_version = as.character(packageVersion("glmnet")))

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("Rescale_FCIR_I_data_N", n, "_P", p, ".Rdata"))
    est_filename = file.path(rdata_dir, paste0("Rescale_FCIR_I_estimates_lasso_sparse_Beta_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("FCIR_I sparse-Beta data not found, skipping:", data_filename))
      next
    }
    if (file.exists(est_filename)) {
      print(paste("FCIR_I sparse-Beta lasso estimates already exist for N =", n, ", P =", p, "- Skipping."))
      next
    }

    # Load into a separate environment: the data file also stores B_reps, N, P,
    # which must not overwrite the settings above.
    d = new.env()
    load(data_filename, envir = d)   # Y, X, Tr, Beta_mat, alpha_0, A_mat, ...
    X = d$X; Tr = d$Tr
    Beta_mat = d$Beta_mat; alpha_0 = d$alpha_0; A_mat = d$A_mat
    print(paste("Starting", B_reps, "FCIR_I sparse-Beta lasso fits for N =", n, "P =", p,
                "with", n_workers, "workers ..."))

    Y_slices = asplit(d$Y[, , seq_len(B_reps), drop = FALSE], 3)
    rm(d)

    cl = parallel::makeCluster(n_workers)
    registerDoParallel(cl)
    parallel::clusterEvalQ(cl, compiler::enableJIT(0))
    parallel::clusterCall(cl, function(dir) {
      source(file.path(dir, "estimate_FCIR_I.r"))
      source(file.path(dir, "estimate_adaptive_lasso_FCIR_I.r"))
    }, dir = fcir_i_code_dir)

    config_start = Sys.time()
    tryCatch({
      fits = foreach(b = seq_along(Y_slices), Y_b = Y_slices,
                     .packages = "glmnet") %dopar% {
        for (attempt in 1:3) {
          out = tryCatch({
            set.seed(seed + b)   # reproducible CV folds
            t0 = proc.time()[["elapsed"]]
            res = estimate_adaptive_lasso_FCIR_I(
              Y = Y_b, X = X, Tr = Tr,
              lambda = lambda_rule, use_cv = TRUE, cv_group_by_site = TRUE,
              adaptive = FALSE)
            list(Beta_mat = res$Beta_mat, alpha_0 = res$alpha_0, A_mat = res$A_mat,
                 lambda = res$selected_lambda,
                 time = proc.time()[["elapsed"]] - t0, attempts = attempt)
          }, error = function(e) if (attempt == 3) stop(e) else NULL)
          if (!is.null(out)) break
        }
        out
      }
    }, finally = {
      try(parallel::stopCluster(cl), silent = TRUE)
      registerDoSEQ()
    })
    wall_time_sec = as.numeric(difftime(Sys.time(), config_start, units = "secs"))
    rm(Y_slices)

    est_Beta_mat = simplify2array(lapply(fits, function(f) f$Beta_mat))
    est_alpha_0 = t(vapply(fits, function(f) f$alpha_0, numeric(L)))
    est_A_mat = simplify2array(lapply(fits, function(f) f$A_mat))
    est_lambda = vapply(fits, function(f) f$lambda, numeric(1))
    est_time_sec = vapply(fits, function(f) f$time, numeric(1))
    est_attempts = vapply(fits, function(f) f$attempts, numeric(1))
    rm(fits)

    penalty = "lasso"
    alpha_value = 1
    cv_type = "grouped_by_site"
    beta_structure = "sparse"
    save(est_Beta_mat, est_alpha_0, est_A_mat, est_lambda, est_time_sec, est_attempts,
         Beta_mat, alpha_0, A_mat,
         N = n, P = p, L, K, B_reps, seed,
         penalty, alpha_value, lambda_rule, cv_type, beta_structure,
         wall_time_sec, n_workers, machine_info,
         file = est_filename)
    print(sprintf("Saved: %s (wall %.1f s, mean %.3f s per fit, %d retried)",
                  est_filename, wall_time_sec, mean(est_time_sec), sum(est_attempts > 1)))
  }
}

print(paste("Total FCIR_I sparse-Beta lasso execution time:", format(Sys.time() - total_start)))
