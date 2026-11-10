# RMSE comparison of the adaptive lasso against the comparison estimators
# (unpenalized MPLE and ordinary lasso for all three models; the FCIR_I and
# FCIR_M lasso files come from main_lasso_sparse_Beta_FCIR_I.r and
# main_lasso_FCIR_M.r and are skipped if they have not been run).
# Also writes the support-recovery metrics of the ordinary lasso for B and A.
#
# RMSE for a parameter block in replicate r is sqrt(mean((est - true)^2)),
# exactly as in the evaluate_selection_* scripts; Theta_int uses its upper
# triangle only. Besides the Monte Carlo mean we keep the median and the number
# of replicates whose estimate diverged (max |est| > diverge_cut), because the
# unpenalized MPLE occasionally runs off to very large values under
# (quasi-)separation and then dominates the mean.

project_root = "D:/ising model thesis/-Fourth-Corner-Ising-Model-"
sim_dir = file.path(project_root, "Simulation_Results")

tol = 1e-8
diverge_cut = 10

# Replicates with non-finite estimates (aliased glm coefficients) are dropped
# and counted in n_nonfinite.
rmse_summary <- function(true_vector, est_matrix) {
  ok <- apply(is.finite(est_matrix), 1, all)
  est_matrix <- est_matrix[ok, , drop = FALSE]
  rmse <- sqrt(rowMeans(sweep(est_matrix, 2, true_vector)^2))
  data.frame(Mean_RMSE = mean(rmse), Median_RMSE = median(rmse), SD_RMSE = sd(rmse),
             n_diverged = sum(apply(abs(est_matrix), 1, max) > diverge_cut),
             n_nonfinite = sum(!ok))
}

selection_summary <- function(true_vector, est_matrix) {
  truth <- abs(true_vector) > tol
  est <- abs(est_matrix) > tol
  tp <- rowSums(sweep(est, 2, truth, "&"))
  tn <- rowSums(sweep(!est, 2, !truth, "&"))
  fp <- rowSums(sweep(est, 2, !truth, "&"))
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall <- tp / sum(truth)
  f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  data.frame(correct_zeros = mean(tn), correct_nonzeros = mean(tp), F1 = mean(f1))
}

flatten_arr <- function(est_arr) t(apply(est_arr, 3, as.vector))
flatten_upper <- function(est_arr) {
  ut <- upper.tri(matrix(0, dim(est_arr)[1], dim(est_arr)[2]))
  t(apply(est_arr, 3, function(m) m[ut]))
}

# blocks: named list of function(env) -> list(true = vector, est = B_reps x d matrix)
evaluate_files <- function(files, blocks, model) {
  out <- list()
  for (i in seq_len(nrow(files))) {
    if (!file.exists(files$file[i])) {
      message("Missing, skipping: ", files$file[i])
      next
    }
    e <- new.env()
    load(files$file[i], envir = e)
    for (b in names(blocks)) {
      x <- blocks[[b]](e)
      out[[length(out) + 1]] <- cbind(
        data.frame(Model = model, Estimator = files$Estimator[i], Block = b,
                   N = files$N[i], P = files$P[i]),
        rmse_summary(x$true, x$est))
    }
  }
  do.call(rbind, out)
}

grid_files <- function(Ns, Ps, patterns) {
  g <- expand.grid(N = Ns, P = Ps, Estimator = names(patterns), stringsAsFactors = FALSE)
  g$file <- mapply(function(n, p, est) sprintf(patterns[[est]], n, p), g$N, g$P, g$Estimator)
  g
}

# ---- Full FCIR ----------------------------------------------------------------
fcir_dir <- file.path(sim_dir, "FCIR", "Rdata_Sparse")
fcir_files <- grid_files(c(50, 100, 200, 400, 800), c(30, 60), list(
  adaptive_lasso = file.path(fcir_dir, "FCIR_estimates_adaptive_lasso_unpenalized_N%d_P%d.Rdata"),
  lasso = file.path(fcir_dir, "FCIR_estimates_lasso_N%d_P%d.Rdata"),
  unpenalized = file.path(fcir_dir, "FCIR_estimates_unpenalized_N%d_P%d.Rdata")))
fcir_blocks <- list(
  beta_0 = function(e) list(true = e$beta_0, est = e$est_beta_0),
  B = function(e) list(true = as.vector(e$B_mat), est = flatten_arr(e$est_B_mat)),
  alpha_0 = function(e) list(true = e$alpha_0, est = e$est_alpha_0),
  A = function(e) list(true = as.vector(e$A_mat), est = flatten_arr(e$est_A_mat)))
fcir_rmse <- evaluate_files(fcir_files, fcir_blocks, "FCIR")

# Ordinary-lasso selection metrics for the sparse fourth-corner matrices.
fcir_lasso_sel <- list()
for (i in which(fcir_files$Estimator == "lasso" & file.exists(fcir_files$file))) {
  e <- new.env()
  load(fcir_files$file[i], envir = e)
  for (b in c("B", "A")) {
    x <- fcir_blocks[[b]](e)
    fcir_lasso_sel[[length(fcir_lasso_sel) + 1]] <- cbind(
      data.frame(Block = b, N = fcir_files$N[i], P = fcir_files$P[i],
                 n_z = sum(abs(x$true) <= tol), n_nz = sum(abs(x$true) > tol)),
      selection_summary(x$true, x$est), rmse_summary(x$true, x$est)["Mean_RMSE"])
  }
}
fcir_lasso_sel <- do.call(rbind, fcir_lasso_sel)

# ---- FCIR_I (sparse Beta, rescaled data used in the thesis) -----------------
fcir_i_dir <- file.path(sim_dir, "FCIR_I", "Rdata_Sparse_Beta")
fcir_i_files <- grid_files(c(50, 100, 200, 400, 800), c(30, 60), list(
  adaptive_lasso = file.path(fcir_i_dir, "Rescale_FCIR_I_estimates_adaptive_lasso_sparse_Beta_unpenalized_N%d_P%d.Rdata"),
  lasso = file.path(fcir_i_dir, "Rescale_FCIR_I_estimates_lasso_sparse_Beta_N%d_P%d.Rdata"),
  unpenalized = file.path(fcir_i_dir, "Rescale_FCIR_I_estimates_unpenalized_sparse_Beta_N%d_P%d.Rdata")))
fcir_i_blocks <- list(
  Beta = function(e) list(true = as.vector(e$Beta_mat), est = flatten_arr(e$est_Beta_mat)),
  alpha_0 = function(e) list(true = e$alpha_0, est = e$est_alpha_0),
  A = function(e) list(true = as.vector(e$A_mat), est = flatten_arr(e$est_A_mat)))
fcir_i_rmse <- evaluate_files(fcir_i_files, fcir_i_blocks, "FCIR_I")

# ---- FCIR_M (L = 3) ------------------------------------------------------------
# The adaptive-lasso estimates behind the thesis tables were not kept, so the
# adaptive-lasso mean RMSE is taken from New_FCIR_M_adaptive_lasso_selection_unpenalized.csv.
fcir_m_dir <- file.path(sim_dir, "FCIR_M", "Rdata")
fcir_m_files <- grid_files(c(50, 100, 200, 400, 800), c(10, 20), list(
  lasso = file.path(fcir_m_dir, "FCIR_M_estimates_lasso_L3_N%d_P%d.Rdata"),
  unpenalized = file.path(fcir_m_dir, "FCIR_M_estimates_unpenalized_L3_N%d_P%d.Rdata")))
fcir_m_blocks <- list(
  beta_0 = function(e) list(true = e$beta_0, est = e$est_beta_0),
  B = function(e) list(true = as.vector(e$B_mat), est = flatten_arr(e$est_B_mat)),
  Theta = function(e) list(true = e$Theta_int[upper.tri(e$Theta_int)], est = flatten_upper(e$est_Theta_int)))
fcir_m_rmse <- evaluate_files(fcir_m_files, fcir_m_blocks, "FCIR_M")

fcir_m_ad <- read.csv(file.path(sim_dir, "FCIR_M", "tables", "New_FCIR_M_adaptive_lasso_selection_unpenalized.csv"),
                      check.names = FALSE)
fcir_m_ad <- fcir_m_ad[fcir_m_ad$Metric == "RMSE", ]
fcir_m_ad <- data.frame(Model = "FCIR_M", Estimator = "adaptive_lasso",
                        Block = sub("Theta_int", "Theta", fcir_m_ad$Matrix),
                        N = fcir_m_ad$N, P = fcir_m_ad$P,
                        Mean_RMSE = fcir_m_ad$Mean, Median_RMSE = NA,
                        SD_RMSE = fcir_m_ad$SD, n_diverged = NA, n_nonfinite = NA)
fcir_m_rmse <- rbind(fcir_m_ad, fcir_m_rmse)

# ---- Write ----------------------------------------------------------------------
write_out <- function(df, model_dir, name) {
  if (is.null(df)) return(invisible())
  est <- if (is.null(df$Estimator)) rep("", nrow(df)) else df$Estimator
  df <- df[order(df$Block, df$P, df$N, est), ]
  f <- file.path(sim_dir, model_dir, "tables", name)
  write.csv(df, f, row.names = FALSE)
  message("Saved: ", f)
  print(df, digits = 4, row.names = FALSE)
}
write_out(fcir_rmse, "FCIR", "estimator_RMSE_comparison.csv")
write_out(fcir_lasso_sel, "FCIR", "lasso_selection_wide.csv")
write_out(fcir_i_rmse, "FCIR_I", "estimator_RMSE_comparison.csv")
write_out(fcir_m_rmse, "FCIR_M", "estimator_RMSE_comparison.csv")
