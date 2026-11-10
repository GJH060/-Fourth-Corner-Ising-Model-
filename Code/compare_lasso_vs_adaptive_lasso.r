# Ordinary lasso vs adaptive lasso, side by side, for all three models:
#   FCIR   : B, A           (+ RMSE of beta_0, alpha_0)
#   FCIR_I : Beta, A        (+ RMSE of alpha_0)   sparse-Beta, rescaled data
#   FCIR_M : B, Theta_int   (+ RMSE of beta_0)    L = 3, B penalised
#
# Metrics per replicate, then averaged over replicates, exactly as in the
# evaluate_selection_* scripts: number of correct zeros, number of correct
# non-zeros, F1 of the support, and RMSE = sqrt(mean((est - true)^2)).
# Theta_int uses its upper triangle only. Beta includes the (unpenalized)
# species-specific intercepts, as in the thesis table.
#
# Inputs: the adaptive lasso estimates of the thesis, and the ordinary lasso
# estimates from FCIR/main_lasso_and_unpenalized.r,
# FCIR_I/main_lasso_sparse_Beta_FCIR_I.r and FCIR_M/main_lasso_FCIR_M.r.
# Note that for FCIR_M the adaptive lasso file is the L3 rerun
# (New_FCIR_M_estimates_adaptive_lasso_L3_*), which reproduces the thesis
# table up to Monte Carlo error in the CV folds.
#
# Outputs (Simulation_Results/lasso_vs_adaptive_lasso/):
#   lasso_vs_adaptive_lasso_long.csv  one row per model/block/N/P/estimator
#   lasso_vs_adaptive_lasso_wide.csv  lasso and adaptive lasso in adjacent columns
#   lasso_vs_adaptive_lasso_wide.txt  the same, as a bordered ASCII table
#   lasso_computation_time.csv        recorded run times of the new lasso fits

project_root = Sys.getenv("FCIR_PROJECT_ROOT",
                          unset = "D:/ising model thesis/-Fourth-Corner-Ising-Model-")
sim_dir = file.path(project_root, "Simulation_Results")
out_dir = file.path(sim_dir, "lasso_vs_adaptive_lasso")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

tol = 1e-8

block_metrics = function(true_vec, est_mat) {
  truth = abs(true_vec) > tol
  est = abs(est_mat) > tol
  tp = rowSums(sweep(est, 2, truth, "&"))
  tn = rowSums(sweep(!est, 2, !truth, "&"))
  fp = rowSums(sweep(est, 2, !truth, "&"))
  fn = sum(truth) - tp
  precision = ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall = ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1 = ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  rmse = sqrt(rowMeans(sweep(est_mat, 2, true_vec)^2))
  data.frame(n_z = sum(!truth), n_nz = sum(truth),
             correct_zeros = mean(tn), correct_nonzeros = mean(tp),
             F1 = mean(f1), Mean_RMSE = mean(rmse), SD_RMSE = sd(rmse))
}

flatten_arr = function(a) t(apply(a, 3, as.vector))
flatten_upper = function(a) {
  ut = upper.tri(matrix(0, dim(a)[1], dim(a)[2]))
  t(apply(a, 3, function(m) m[ut]))
}

specs = list(
  FCIR = list(
    Ns = c(50, 100, 200, 400, 800), Ps = c(30, 60),
    files = list(
      adaptive_lasso = file.path(sim_dir, "FCIR", "Rdata_Sparse", "FCIR_estimates_adaptive_lasso_unpenalized_N%d_P%d.Rdata"),
      lasso = file.path(sim_dir, "FCIR", "Rdata_Sparse", "FCIR_estimates_lasso_N%d_P%d.Rdata")),
    blocks = list(
      B = function(e) list(true = as.vector(e$B_mat), est = flatten_arr(e$est_B_mat)),
      A = function(e) list(true = as.vector(e$A_mat), est = flatten_arr(e$est_A_mat)),
      beta_0 = function(e) list(true = e$beta_0, est = e$est_beta_0),
      alpha_0 = function(e) list(true = e$alpha_0, est = e$est_alpha_0))),
  FCIR_I = list(
    Ns = c(50, 100, 200, 400, 800), Ps = c(30, 60),
    files = list(
      adaptive_lasso = file.path(sim_dir, "FCIR_I", "Rdata_Sparse_Beta", "Rescale_FCIR_I_estimates_adaptive_lasso_sparse_Beta_unpenalized_N%d_P%d.Rdata"),
      lasso = file.path(sim_dir, "FCIR_I", "Rdata_Sparse_Beta", "Rescale_FCIR_I_estimates_lasso_sparse_Beta_N%d_P%d.Rdata")),
    blocks = list(
      Beta = function(e) list(true = as.vector(e$Beta_mat), est = flatten_arr(e$est_Beta_mat)),
      A = function(e) list(true = as.vector(e$A_mat), est = flatten_arr(e$est_A_mat)),
      alpha_0 = function(e) list(true = e$alpha_0, est = e$est_alpha_0))),
  FCIR_M = list(
    Ns = c(50, 100, 200, 400, 800), Ps = c(10, 20),
    files = list(
      adaptive_lasso = file.path(sim_dir, "FCIR_M", "Rdata", "New_FCIR_M_estimates_adaptive_lasso_L3_unpenalized_N%d_P%d.Rdata"),
      lasso = file.path(sim_dir, "FCIR_M", "Rdata", "FCIR_M_estimates_lasso_L3_N%d_P%d.Rdata")),
    blocks = list(
      B = function(e) list(true = as.vector(e$B_mat), est = flatten_arr(e$est_B_mat)),
      Theta_int = function(e) list(true = e$Theta_int[upper.tri(e$Theta_int)], est = flatten_upper(e$est_Theta_int)),
      beta_0 = function(e) list(true = e$beta_0, est = e$est_beta_0)))
)

rows = list()
times = list()
for (m in names(specs)) {
  s = specs[[m]]
  for (n in s$Ns) for (p in s$Ps) for (est in names(s$files)) {
    f = sprintf(s$files[[est]], n, p)
    if (!file.exists(f)) { message("Missing, skipping: ", f); next }
    e = new.env(); load(f, envir = e)
    for (b in names(s$blocks)) {
      x = s$blocks[[b]](e)
      rows[[length(rows) + 1]] = cbind(data.frame(Model = m, Block = b, N = n, P = p, Estimator = est),
                                       block_metrics(x$true, x$est))
    }
    if (est == "lasso" && !is.null(e$est_time_sec)) {
      times[[length(times) + 1]] = data.frame(
        Model = m, N = n, P = p, n_reps = length(e$est_time_sec),
        mean_sec_per_fit = mean(e$est_time_sec), median_sec_per_fit = median(e$est_time_sec),
        wall_time_min = e$wall_time_sec / 60, n_workers = e$n_workers,
        cpu = e$machine_info$cpu)
    }
  }
}
long = do.call(rbind, rows)
long = long[order(match(long$Model, names(specs)), long$Block, long$P, long$N, long$Estimator), ]

key = c("Model", "Block", "N", "P", "n_z", "n_nz")
vals = c("correct_zeros", "correct_nonzeros", "F1", "Mean_RMSE")
wide = merge(long[long$Estimator == "lasso", c(key, vals)],
             long[long$Estimator == "adaptive_lasso", c(key, vals)],
             by = key, suffixes = c("_lasso", "_adaptive"), all = TRUE)
wide = wide[order(match(wide$Model, names(specs)), wide$Block, wide$P, wide$N),
            c(key, as.vector(t(outer(vals, c("_lasso", "_adaptive"), paste0))))]
# Selection metrics are meaningless for the intercept-type blocks (all non-zero
# or unpenalized); keep only their RMSE.
dense = wide$Block %in% c("beta_0", "alpha_0")
wide[dense, grep("^(correct|F1)", names(wide))] = NA

format_ascii_table = function(x, digits = 3) {
  f = as.data.frame(lapply(x, function(col) {
    out = if (is.numeric(col)) format(round(col, digits), trim = TRUE, scientific = FALSE) else as.character(col)
    out[is.na(col)] = ""
    out
  }), stringsAsFactors = FALSE, check.names = FALSE)
  w = pmax(nchar(names(f)), vapply(f, function(col) max(nchar(col)), integer(1)))
  border = paste0("+", paste(strrep("-", w + 2), collapse = "+"), "+")
  line = function(v) paste0("| ", paste(sprintf(paste0("%-", w, "s"), v), collapse = " | "), " |")
  c(border, line(names(f)), border, apply(f, 1, line), border)
}

write.csv(long, file.path(out_dir, "lasso_vs_adaptive_lasso_long.csv"), row.names = FALSE)
write.csv(wide, file.path(out_dir, "lasso_vs_adaptive_lasso_wide.csv"), row.names = FALSE)
writeLines(format_ascii_table(wide), file.path(out_dir, "lasso_vs_adaptive_lasso_wide.txt"))
writeLines(format_ascii_table(wide))
if (length(times) > 0) {
  times = do.call(rbind, times)
  write.csv(times, file.path(out_dir, "lasso_computation_time.csv"), row.names = FALSE)
  print(times, digits = 3, row.names = FALSE)
}
message("Saved tables to ", out_dir)
