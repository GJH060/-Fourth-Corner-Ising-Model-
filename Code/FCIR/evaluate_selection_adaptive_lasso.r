
library(dplyr)

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"

use_dense = FALSE        # FALSE -> Rdata_Sparse, TRUE -> Rdata_Dense
init_method = "unpenalized"  # must match the init used in main_adaptive_lasso.r ("ridge" / "unpenalized")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR",
                      if (use_dense) "Rdata_Dense" else "Rdata_Sparse")
infix = if (use_dense) "_Dense" else ""
out_dir = file.path(project_root, "Simulation_Results", "FCIR", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

tol = 1e-8   # |coef| <= tol is treated as a (selected-out) zero

# true_vector : length-M vector of true coefficients (some exactly zero)
# est_matrix  : B_reps x M matrix, one estimated vector per replicate
evaluate_selection <- function(true_vector, est_matrix, tol = 1e-8) {
  true_support <- as.numeric(abs(true_vector) > tol)
  est_support  <- matrix(as.numeric(abs(est_matrix) > tol),
                         nrow = nrow(est_matrix), ncol = ncol(est_matrix))

  num_sim <- nrow(est_matrix)
  n_correct_zeros <- n_correct_nonzeros <- f1 <- rmse <- numeric(num_sim)

  for (i in 1:num_sim) {
    y_true <- true_support
    y_pred <- est_support[i, ]

    tp <- sum(y_true == 1 & y_pred == 1)
    tn <- sum(y_true == 0 & y_pred == 0)
    fp <- sum(y_true == 0 & y_pred == 1)
    fn <- sum(y_true == 1 & y_pred == 0)

    # Absolute selection counts for zero and non-zero coefficients.
    n_correct_zeros[i]    <- tn
    n_correct_nonzeros[i] <- tp

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1[i]     <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0

    rmse[i]   <- sqrt(mean((est_matrix[i, ] - true_vector)^2))
  }

  data.frame(
    Metric = c("no of correct zeros", "no of correct non-zeros", "F1", "RMSE"),
    Mean   = c(mean(n_correct_zeros), mean(n_correct_nonzeros),
               mean(f1, na.rm = TRUE), mean(rmse, na.rm = TRUE)),
    SD     = c(sd(n_correct_zeros), sd(n_correct_nonzeros),
               sd(f1, na.rm = TRUE), sd(rmse, na.rm = TRUE))
  )
}

# Dense blocks have no true zeros, so FPR/FDR/F1 are not meaningful here.
# nonzero_recovery_rate is the proportion of truly nonzero coefficients that
# remain nonzero after adaptive-lasso shrinkage.
evaluate_dense_block <- function(true_vector, est_matrix, tol = 1e-8) {
  true_nonzero <- abs(true_vector) > tol
  num_sim <- nrow(est_matrix)
  nonzero_recovery_rate <- rmse <- numeric(num_sim)

  for (i in 1:num_sim) {
    est_nonzero <- abs(est_matrix[i, ]) > tol
    nonzero_recovery_rate[i] <- mean(est_nonzero[true_nonzero])
    rmse[i] <- sqrt(mean((est_matrix[i, ] - true_vector)^2))
  }

  data.frame(
    Metric = c("nonzero_recovery_rate", "RMSE"),
    Mean = c(mean(nonzero_recovery_rate, na.rm = TRUE), mean(rmse, na.rm = TRUE)),
    SD = c(sd(nonzero_recovery_rate, na.rm = TRUE), sd(rmse, na.rm = TRUE))
  )
}

# Flatten an L x K x B_reps array into a B_reps x (L*K) matrix (column-major per rep,
# matching as.vector of an L x K matrix so columns line up with the true vector).
flatten_estimates <- function(est_arr) {
  B_reps <- dim(est_arr)[3]
  t(apply(est_arr, 3, as.vector))   # B_reps x (L*K)
}

format_ascii_table <- function(x, digits = 3) {
  x <- as.data.frame(x)
  formatted <- lapply(x, function(col) {
    if (is.numeric(col)) {
      out <- format(round(col, digits), trim = TRUE, scientific = FALSE)
    } else {
      out <- as.character(col)
    }
    out[is.na(out)] <- ""
    out
  })
  formatted <- as.data.frame(formatted, stringsAsFactors = FALSE, check.names = FALSE)
  names(formatted) <- names(x)

  widths <- pmax(
    nchar(names(formatted)),
    vapply(formatted, function(col) max(nchar(col)), integer(1))
  )
  border <- paste0("+", paste(vapply(widths, function(w) paste(rep("-", w + 2), collapse = ""), character(1)), collapse = "+"), "+")
  row_line <- function(vals) {
    paste0("| ", paste(mapply(function(val, w) sprintf(paste0("%-", w, "s"), val), vals, widths), collapse = " | "), " |")
  }

  c(
    border,
    row_line(names(formatted)),
    border,
    apply(formatted, 1, row_line),
    border
  )
}

results = list()

for (n in Ns) {
  for (p in Ps) {
    f = file.path(rdata_dir, paste0("FCIR_estimates_adaptive_lasso_", init_method, infix, "_N", n, "_P", p, ".Rdata"))
    if (!file.exists(f)) {
      warning(paste("Missing:", f, "- run main_adaptive_lasso.r first. Skipping."))
      next
    }
    e = new.env(); load(f, envir = e)

    for (mat in c("B", "A")) {
      true_vec = as.vector(if (mat == "B") e$B_mat else e$A_mat)
      est_mat  = flatten_estimates(if (mat == "B") e$est_B_mat else e$est_A_mat)
      summ = evaluate_selection(true_vec, est_mat, tol = tol)
      summ$N = n; summ$P = p; summ$Matrix = mat
      summ$n_nonzero = sum(abs(true_vec) > tol)
      summ$n_total   = length(true_vec)
      results[[paste(mat, n, p, sep = "_")]] = summ
    }

    for (block in c("beta_0", "alpha_0")) {
      true_vec = as.vector(if (block == "beta_0") e$beta_0 else e$alpha_0)
      est_mat = if (block == "beta_0") e$est_beta_0 else e$est_alpha_0
      summ = evaluate_dense_block(true_vec, est_mat, tol = tol)
      summ$N = n; summ$P = p; summ$Matrix = block
      summ$n_nonzero = sum(abs(true_vec) > tol)
      summ$n_total = length(true_vec)
      results[[paste(block, n, p, sep = "_")]] = summ
    }
    message("Evaluated adaptive lasso metrics: N=", n, ", P=", p)
  }
}

if (length(results) == 0) stop("No adaptive lasso estimate files found.")

df = bind_rows(results) %>%
  select(Matrix, N, P, n_nonzero, n_total, Metric, Mean, SD) %>%
  arrange(Matrix, P, N, Metric)

scenario_tag = if (use_dense) "Dense" else "Sparse"
out_csv = file.path(out_dir, paste0("adaptive_lasso_selection_", init_method, "_", scenario_tag, ".csv"))
write.csv(df, out_csv, row.names = FALSE)
message("Saved selection summary table: ", out_csv)

# Wide selection table for sparse fourth-corner matrices only.
selection_wide = df %>%
  filter(Matrix %in% c("B", "A")) %>%
  select(-SD) %>%
  tidyr::pivot_wider(id_cols = c(Matrix, N, P, n_nonzero, n_total),
                     names_from = Metric, values_from = Mean) %>%
  select(Matrix, N, P, n_nonzero, n_total,
         `no of correct zeros`, `no of correct non-zeros`, F1, RMSE) %>%
  arrange(Matrix, P, N)
selection_wide_csv = file.path(out_dir, paste0("adaptive_lasso_selection_wide_", init_method, "_", scenario_tag, ".csv"))
write.csv(selection_wide, selection_wide_csv, row.names = FALSE)
selection_wide_txt = file.path(out_dir, paste0("adaptive_lasso_selection_wide_", init_method, "_", scenario_tag, ".txt"))
writeLines(format_ascii_table(selection_wide), selection_wide_txt)

# Wide diagnostic table for dense baseline blocks only.
dense_wide = df %>%
  filter(Matrix %in% c("beta_0", "alpha_0")) %>%
  select(-SD) %>%
  tidyr::pivot_wider(id_cols = c(Matrix, N, P, n_nonzero, n_total),
                     names_from = Metric, values_from = Mean) %>%
  select(Matrix, N, P, n_nonzero, n_total, nonzero_recovery_rate, RMSE) %>%
  arrange(Matrix, P, N)
dense_wide_csv = file.path(out_dir, paste0("adaptive_lasso_dense_diagnostics_wide_", init_method, "_", scenario_tag, ".csv"))
write.csv(dense_wide, dense_wide_csv, row.names = FALSE)
dense_wide_txt = file.path(out_dir, paste0("adaptive_lasso_dense_diagnostics_wide_", init_method, "_", scenario_tag, ".txt"))
writeLines(format_ascii_table(dense_wide), dense_wide_txt)

message("Sparse B/A selection table:")
writeLines(format_ascii_table(selection_wide))
message("Dense beta_0/alpha_0 diagnostic table:")
writeLines(format_ascii_table(dense_wide))
message("Saved wide selection summary table: ", selection_wide_csv)
message("Saved bordered wide selection table: ", selection_wide_txt)
message("Saved dense diagnostic summary table: ", dense_wide_csv)
message("Saved bordered dense diagnostic table: ", dense_wide_txt)
message("Finished adaptive lasso selection evaluation.")

