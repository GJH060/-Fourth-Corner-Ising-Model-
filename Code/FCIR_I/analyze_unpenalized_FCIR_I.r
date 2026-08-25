library(ggplot2)
library(dplyr)

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"

# Must match the setting used in main_unpenalized_FCIR_I.r.
use_sparse_beta = TRUE

# Estimates fitted on the rescaled data carry a "Rescale_" prefix; the same
# prefix goes on the tables and plots so the earlier non-rescaled output is
# left intact. Set to "" to analyse the original estimates instead.
est_prefix = "Rescale_"

if (use_sparse_beta) {
  rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata_Sparse_Beta")
  est_tag = "unpenalized_sparse_Beta"
  plot_subdir = "unpenalized_sparse_Beta"
  data_label = "FCIR_I sparse-Beta"
  title_label = "FCIR_I Sparse-Beta Unpenalized"
} else {
  rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata")
  est_tag = "unpenalized"
  plot_subdir = "unpenalized"
  data_label = "FCIR_I"
  title_label = "FCIR_I Unpenalized"
}

plot_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "plots",
                     plot_subdir)
table_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "tables")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

flatten_estimates <- function(est_arr) {
  t(apply(est_arr, 3, as.vector))
}

calculate_rmse <- function(est_matrix, true_vector) {
  error_matrix = sweep(as.matrix(est_matrix), 2, as.numeric(true_vector), "-")
  sqrt(rowMeans(error_matrix^2))
}

tol = 1e-8

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
beta_scatter_results = list()

for (n in Ns) {
  for (p in Ps) {
    filename = file.path(
      rdata_dir,
      paste0(est_prefix, "FCIR_I_estimates_", est_tag,
             "_N", n, "_P", p, ".Rdata")
    )

    if (!file.exists(filename)) {
      warning(paste("File not found:", filename, "- Skipping."))
      next
    }

    e = new.env()
    load(filename, envir = e)

    block_inputs = list(
      Beta_mat = list(
        estimates = flatten_estimates(e$est_Beta_mat),
        truth = as.vector(e$Beta_mat)
      ),
      alpha_0 = list(
        estimates = e$est_alpha_0,
        truth = as.vector(e$alpha_0)
      ),
      A_mat = list(
        estimates = flatten_estimates(e$est_A_mat),
        truth = as.vector(e$A_mat)
      )
    )

    beta_rmse = calculate_rmse(
      block_inputs$Beta_mat$estimates,
      block_inputs$Beta_mat$truth
    )
    beta_scatter_results[[paste(n, p, sep = "_")]] = data.frame(
      N = n,
      P = p,
      Replicate = seq_along(beta_rmse),
      Max_Abs_Beta = apply(abs(e$est_Beta_mat), 3, max, na.rm = TRUE),
      RMSE = beta_rmse
    )

    for (block in names(block_inputs)) {
      rmse = calculate_rmse(
        block_inputs[[block]]$estimates,
        block_inputs[[block]]$truth
      )

      truth = block_inputs[[block]]$truth
      results[[paste(block, n, p, sep = "_")]] = data.frame(
        Matrix = block,
        N = n,
        P = p,
        n_nonzero = sum(abs(truth) > tol),
        n_total = length(truth),
        Replicate = seq_along(rmse),
        RMSE = rmse
      )
    }

    message("Evaluated ", data_label, " unpenalized RMSE: N=", n, ", P=", p)
  }
}

if (length(results) == 0) {
  stop(paste("No", data_label, "unpenalized estimate files were found."))
}

rmse_df = bind_rows(results)

summary_df = rmse_df %>%
  group_by(Matrix, N, P, n_nonzero, n_total) %>%
  summarise(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    SD_RMSE = sd(RMSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Matrix, P, N)

summary_csv = file.path(
  table_dir,
  paste0(est_prefix, "FCIR_I_", est_tag, "_RMSE_summary.csv")
)
write.csv(summary_df, summary_csv, row.names = FALSE)

rmse_wide = summary_df %>%
  select(Matrix, N, P, n_nonzero, n_total, Mean_RMSE, SD_RMSE) %>%
  arrange(Matrix, P, N)

rmse_wide_csv = file.path(
  table_dir,
  paste0(est_prefix, "FCIR_I_", est_tag, "_RMSE_wide.csv")
)
write.csv(rmse_wide, rmse_wide_csv, row.names = FALSE)

rmse_wide_txt = file.path(
  table_dir,
  paste0(est_prefix, "FCIR_I_", est_tag, "_RMSE_wide.txt")
)
writeLines(format_ascii_table(rmse_wide), rmse_wide_txt)

message(data_label, " unpenalized RMSE table:")
writeLines(format_ascii_table(rmse_wide))

rmse_df$N = factor(rmse_df$N, levels = Ns)
rmse_df$P = factor(rmse_df$P)

rmse_plot_df = rmse_df %>%
  filter(is.finite(RMSE), RMSE > 0)

plot_obj = ggplot(rmse_plot_df, aes(x = N, y = RMSE)) +
  geom_boxplot(fill = "skyblue", color = "black",
               outlier.shape = 16, outlier.alpha = 0.3, width = 0.55) +
  scale_y_log10() +
  facet_grid(Matrix ~ P, scales = "free_y",
             labeller = labeller(P = function(x) paste0("P = ", x))) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste(title_label, "RMSE Convergence (Log Scale)"),
    x = "Sample Size (N)",
    y = "RMSE (log10 scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold")
  )

plot_file = file.path(
  plot_dir,
  paste0(est_prefix, "FCIR_I_", est_tag, "_RMSE_convergence.png")
)
ggsave(plot_file, plot_obj, width = 11, height = 9, dpi = 300)

beta_scatter_df = bind_rows(beta_scatter_results) %>%
  filter(is.finite(Max_Abs_Beta), Max_Abs_Beta > 0,
         is.finite(RMSE), RMSE > 0) %>%
  mutate(
    N = factor(N, levels = Ns),
    P = factor(P)
  )

beta_scatter_plot = ggplot(
  beta_scatter_df,
  aes(x = Max_Abs_Beta, y = RMSE)
) +
  geom_point(alpha = 0.35, size = 1.2, color = "steelblue") +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(N ~ P, labeller = labeller(
    N = function(x) paste0("N = ", x),
    P = function(x) paste0("P = ", x)
  )) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Maximum Absolute Beta Estimate vs Beta RMSE",
    subtitle = "Each point represents one Monte Carlo replicate",
    x = "max |estimated Beta| (log10 scale)",
    y = "Beta RMSE (log10 scale)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold")
  )

beta_scatter_file = file.path(
  plot_dir,
  paste0(est_prefix, "FCIR_I_", est_tag, "_max_abs_Beta_vs_RMSE.png")
)
ggsave(beta_scatter_file, beta_scatter_plot,
       width = 9, height = 12, dpi = 300)

message("Saved RMSE summary: ", summary_csv)
message("Saved wide RMSE summary table: ", rmse_wide_csv)
message("Saved bordered wide RMSE table: ", rmse_wide_txt)
message("Saved log-scale RMSE convergence plot: ", plot_file)
message("Saved max |Beta| versus RMSE scatter plot: ", beta_scatter_file)
message("Finished ", data_label, " unpenalized RMSE analysis.")
