# Retrospective computation-time benchmark for the Chapter 4 simulations.
#
# The original simulation runs did not record computation time. This script
# re-times the estimators on the same simulated data with the same code and
# settings, for a subset of replicates per configuration:
#
#   models     : FCIR (Rdata_Sparse), FCIR_I (sparse Beta, rescaled), FCIR_M (L = 3)
#   estimators : unpenalized MPLE, ordinary lasso, adaptive lasso
#                (the adaptive lasso time includes its unpenalized initial fit)
#
# Replicates 1..n_bench of each configuration are fitted in parallel, one fit
# per worker, the same way the simulation scripts ran them. Each fit is timed
# inside its worker (elapsed seconds), so the per-fit times reflect a fully
# loaded machine, as in the original runs. The projected time for the full
# 1000 replicates is mean_sec * 1000 / n_workers.
#
# Results are appended to Simulation_Results/computation_time/ one row per
# configuration as soon as it finishes, so an interrupted run resumes where it
# stopped. Delete the CSV to start over.
#
# Do not run this while other simulations are running: they compete for the
# same cores and inflate the times.

project_root = Sys.getenv("FCIR_PROJECT_ROOT",
                          unset = "D:/ising model thesis/-Fourth-Corner-Ising-Model-")
code_dir = file.path(project_root, "Code")
sim_dir = file.path(project_root, "Simulation_Results")
out_dir = file.path(sim_dir, "computation_time")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
fits_csv = file.path(out_dir, "benchmark_fit_times.csv")       # one row per fit
summary_csv = file.path(out_dir, "benchmark_summary.csv")      # one row per configuration

library(foreach)
library(doParallel)

n_workers = max(1, parallel::detectCores() - 1)
n_bench = as.integer(Sys.getenv("BENCH_REPS", unset = n_workers))   # replicates timed per configuration
B_reps_full = 1000
seed = 2026

models = list(
  FCIR = list(
    code = c("FCIR/estimate_FCIR.r", "FCIR/estimate_adaptive_lasso_FCIR.r"),
    data = function(n, p) file.path(sim_dir, "FCIR", "Rdata_Sparse", sprintf("FCIR_data_N%d_P%d.Rdata", n, p)),
    Ns = c(50, 100, 200, 400, 800), Ps = c(30, 60),
    fit = list(
      unpenalized = function(Y, X, Tr) estimate_unpenalized_FCIR(Y = Y, X = X, Tr = Tr, standardize = FALSE),
      lasso = function(Y, X, Tr) estimate_adaptive_lasso_FCIR(Y = Y, X = X, Tr = Tr,
        lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE, adaptive = FALSE),
      adaptive_lasso = function(Y, X, Tr) estimate_adaptive_lasso_FCIR(Y = Y, X = X, Tr = Tr,
        gamma = 1, init = "unpenalized", lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE))),
  FCIR_I = list(
    code = c("FCIR_I/estimate_FCIR_I.r", "FCIR_I/estimate_adaptive_lasso_FCIR_I.r"),
    data = function(n, p) file.path(sim_dir, "FCIR_I", "Rdata_Sparse_Beta", sprintf("Rescale_FCIR_I_data_N%d_P%d.Rdata", n, p)),
    Ns = c(50, 100, 200, 400, 800), Ps = c(30, 60),
    fit = list(
      unpenalized = function(Y, X, Tr) estimate_unpenalized_FCIR_I(Y = Y, X = X, Tr = Tr),
      lasso = function(Y, X, Tr) estimate_adaptive_lasso_FCIR_I(Y = Y, X = X, Tr = Tr,
        lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE, adaptive = FALSE),
      adaptive_lasso = function(Y, X, Tr) estimate_adaptive_lasso_FCIR_I(Y = Y, X = X, Tr = Tr,
        gamma = 1, init = "unpenalized", lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE))),
  FCIR_M = list(
    code = c("FCIR_M/estimate_FCIR_M.r", "FCIR_M/estimate_adaptive_lasso_FCIR_M.r"),
    data = function(n, p) file.path(sim_dir, "FCIR_M", "Rdata", sprintf("FCIR_M_data_L3_N%d_P%d.Rdata", n, p)),
    Ns = c(50, 100, 200, 400, 800), Ps = c(10, 20),
    fit = list(
      unpenalized = function(Y, X, Tr) estimate_unpenalized_FCIR_M(Y = Y, X = X, Tr = Tr),
      lasso = function(Y, X, Tr) estimate_adaptive_lasso_FCIR_M(Y = Y, X = X, Tr = Tr,
        lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE, unpenalize_B = FALSE, adaptive = FALSE),
      adaptive_lasso = function(Y, X, Tr) estimate_adaptive_lasso_FCIR_M(Y = Y, X = X, Tr = Tr,
        gamma = 1, init = "unpenalized", lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE,
        unpenalize_B = FALSE)))
)

# Optional filters, e.g. BENCH_MODELS="FCIR_M" BENCH_ESTIMATORS="lasso,adaptive_lasso"
only = function(var, all) if (nzchar(Sys.getenv(var))) intersect(all, strsplit(Sys.getenv(var), ",")[[1]]) else all
model_names = only("BENCH_MODELS", names(models))
estimator_names = only("BENCH_ESTIMATORS", c("unpenalized", "lasso", "adaptive_lasso"))

done = if (file.exists(summary_csv)) read.csv(summary_csv, stringsAsFactors = FALSE) else NULL
is_done = function(m, est, n, p) {
  !is.null(done) && any(done$Model == m & done$Estimator == est & done$N == n & done$P == p)
}
append_csv = function(df, f) {
  write.table(df, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
}

cpu = Sys.getenv("PROCESSOR_IDENTIFIER")
message("Benchmark: ", n_bench, " replicates per configuration, ", n_workers, " workers, CPU: ", cpu)

for (m in model_names) {
  spec = models[[m]]
  cl = parallel::makeCluster(n_workers)
  registerDoParallel(cl)
  parallel::clusterEvalQ(cl, compiler::enableJIT(0))
  parallel::clusterCall(cl, function(dir, files) {
    for (f in files) source(file.path(dir, f))
  }, dir = code_dir, files = spec$code)

  for (est in estimator_names) {
    fit_fun = spec$fit[[est]]
    for (n in spec$Ns) for (p in spec$Ps) {
      if (is_done(m, est, n, p)) { message("Done already: ", m, " ", est, " N=", n, " P=", p); next }
      f = spec$data(n, p)
      if (!file.exists(f)) { message("Missing data, skipping: ", f); next }

      d = new.env(); load(f, envir = d)
      X = d$X; Tr = d$Tr
      Y_slices = asplit(d$Y[, , seq_len(n_bench), drop = FALSE], 3)
      rm(d)

      batch_start = Sys.time()
      res = tryCatch(
        foreach(b = seq_along(Y_slices), Y_b = Y_slices, .packages = c("glmnet"),
                .export = "fit_fun", .errorhandling = "pass") %dopar% {
          set.seed(seed + b)   # same CV folds as the simulation scripts
          t0 = proc.time()[["elapsed"]]
          suppressWarnings(fit_fun(Y_b, X, Tr))
          proc.time()[["elapsed"]] - t0
        },
        error = function(e) { message("Batch failed: ", conditionMessage(e)); NULL })
      batch_wall = as.numeric(difftime(Sys.time(), batch_start, units = "secs"))
      rm(Y_slices)
      if (is.null(res)) next

      ok = vapply(res, is.numeric, logical(1))
      times = unlist(res[ok])
      append_csv(data.frame(Model = m, Estimator = est, N = n, P = p,
                            replicate = which(ok), elapsed_sec = times), fits_csv)
      row = data.frame(
        Model = m, Estimator = est, N = n, P = p,
        n_timed = length(times), n_failed = sum(!ok),
        mean_sec = mean(times), median_sec = median(times), sd_sec = sd(times),
        min_sec = min(times), max_sec = max(times),
        batch_wall_sec = batch_wall, n_workers = n_workers,
        projected_wall_hours_1000 = mean(times) * B_reps_full / n_workers / 3600,
        projected_cpu_hours_1000 = mean(times) * B_reps_full / 3600,
        cpu = cpu, r_version = R.version.string, date = format(Sys.time(), "%Y-%m-%d %H:%M"))
      append_csv(row, summary_csv)
      message(sprintf("%-7s %-15s N=%-4d P=%-3d mean %.3f s/fit, projected %.2f h for 1000 reps on %d workers",
                      m, est, n, p, row$mean_sec, row$projected_wall_hours_1000, n_workers))
    }
  }
  try(parallel::stopCluster(cl), silent = TRUE)
  registerDoSEQ()
}

# Wide table: mean seconds per fit for each estimator side by side, plus the
# projected wall time of 1000 replicates on n_workers workers.
s = read.csv(summary_csv, stringsAsFactors = FALSE)
s$Estimator = factor(s$Estimator, c("unpenalized", "lasso", "adaptive_lasso"))
s = s[order(match(s$Model, names(models)), s$P, s$N, s$Estimator), ]
wide = reshape(s[, c("Model", "N", "P", "Estimator", "mean_sec", "projected_wall_hours_1000")],
               idvar = c("Model", "N", "P"), timevar = "Estimator", direction = "wide")
names(wide) = sub("^mean_sec\\.", "sec_per_fit_", sub("^projected_wall_hours_1000\\.", "hours_1000reps_", names(wide)))
write.csv(wide, file.path(out_dir, "benchmark_summary_wide.csv"), row.names = FALSE)

message("Saved: ", summary_csv, ", ", fits_csv, " and benchmark_summary_wide.csv")
