library(dplyr)

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata_Sparse_Beta")
out_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

init_method = "unpenalized"
tol = 1e-8

evaluate_selection <- function(true_vector, est_matrix, tol = 1e-8) {
  true_support <- as.numeric(abs(true_vector) > tol)
  est_support <- matrix(as.numeric(abs(est_matrix) > tol),
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

    n_correct_zeros[i] <- tn
    n_correct_nonzeros[i] <- tp

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1[i] <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0

    rmse[i] <- sqrt(mean((est_matrix[i, ] - true_vector)^2))
  }

  data.frame(
    Metric = c("no of correct zeros", "no of correct non-zeros", "F1", "RMSE"),
    Mean = c(mean(n_correct_zeros), mean(n_correct_nonzeros),
             mean(f1, na.rm = TRUE), mean(rmse, na.rm = TRUE)),
    SD = c(sd(n_correct_zeros), sd(n_correct_nonzeros),
           sd(f1, na.rm = TRUE), sd(rmse, na.rm = TRUE))
  )
}

flatten_estimates <- function(est_arr) {
  t(apply(est_arr, 3, as.vector))
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
    f = file.path(
      rdata_dir,
      paste0("FCIR_I_estimates_adaptive_lasso_sparse_Beta_", init_method,
             "_N", n, "_P", p, ".Rdata")
    )
    if (!file.exists(f)) {
      warning(paste("Missing:", f,
                    "- run main_adaptive_lasso_sparse_Beta_FCIR_I.r first. Skipping."))
      next
    }
    e = new.env(); load(f, envir = e)

    true_vec = as.vector(e$Beta_mat)
    est_mat = flatten_estimates(e$est_Beta_mat)
    summ = evaluate_selection(true_vec, est_mat, tol = tol)
    summ$N = n; summ$P = p; summ$Matrix = "Beta_mat"
    summ$n_nonzero = sum(abs(true_vec) > tol)
    summ$n_total = length(true_vec)
    results[[paste("Beta_mat", n, p, sep = "_")]] = summ

    true_vec = as.vector(e$A_mat)
    est_mat = flatten_estimates(e$est_A_mat)
    summ = evaluate_selection(true_vec, est_mat, tol = tol)
    summ$N = n; summ$P = p; summ$Matrix = "A"
    summ$n_nonzero = sum(abs(true_vec) > tol)
    summ$n_total = length(true_vec)
    results[[paste("A", n, p, sep = "_")]] = summ

    true_vec = as.vector(e$alpha_0)
    est_mat = e$est_alpha_0
    summ = evaluate_selection(true_vec, est_mat, tol = tol)
    summ$N = n; summ$P = p; summ$Matrix = "alpha_0"
    summ$n_nonzero = sum(abs(true_vec) > tol)
    summ$n_total = length(true_vec)
    results[[paste("alpha_0", n, p, sep = "_")]] = summ

    message("Evaluated FCIR_I sparse-Beta adaptive lasso metrics: N=", n, ", P=", p)
  }
}

if (length(results) == 0) {
  stop("No FCIR_I sparse-Beta adaptive lasso estimate files found.")
}

df = bind_rows(results) %>%
  select(Matrix, N, P, n_nonzero, n_total, Metric, Mean, SD) %>%
  arrange(Matrix, N, P, Metric)

out_csv = file.path(
  out_dir,
  paste0("FCIR_I_adaptive_lasso_selection_sparse_Beta_", init_method, ".csv")
)
write.csv(df, out_csv, row.names = FALSE)
message("Saved FCIR_I sparse-Beta selection summary table: ", out_csv)

selection_wide = df %>%
  select(-SD) %>%
  tidyr::pivot_wider(id_cols = c(Matrix, N, P, n_nonzero, n_total),
                     names_from = Metric, values_from = Mean) %>%
  select(Matrix, N, P, n_nonzero, n_total,
         `no of correct zeros`, `no of correct non-zeros`, F1, RMSE) %>%
  arrange(Matrix, N, P)

selection_wide_csv = file.path(
  out_dir,
  paste0("FCIR_I_adaptive_lasso_selection_wide_sparse_Beta_", init_method, ".csv")
)
write.csv(selection_wide, selection_wide_csv, row.names = FALSE)

selection_wide_txt = file.path(
  out_dir,
  paste0("FCIR_I_adaptive_lasso_selection_wide_sparse_Beta_", init_method, ".txt")
)
writeLines(format_ascii_table(selection_wide), selection_wide_txt)

message("FCIR_I sparse-Beta adaptive lasso selection table:")
writeLines(format_ascii_table(selection_wide))
message("Saved FCIR_I sparse-Beta wide selection summary table: ", selection_wide_csv)
message("Saved FCIR_I sparse-Beta bordered wide selection table: ", selection_wide_txt)
message("Finished FCIR_I sparse-Beta adaptive lasso selection evaluation.")
