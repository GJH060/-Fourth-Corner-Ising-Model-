# Compare interaction-network selection metrics:
#   Ising  Theta      vs  FCIR_M  Theta_int
# at matching (N, P), for correct zeros / correct nonzeros / F1 / RMSE.

library(ggplot2)
library(dplyr)
library(tidyr)

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"

ising_txt = file.path(
  project_root, "Simulation_Results", "Ising", "tables",
  "Ising_adpative_lasso_selection.txt"
)
fcir_m_csv = file.path(
  project_root, "Simulation_Results", "FCIR_M", "tables",
  "New_FCIR_M_adaptive_lasso_selection_wide_unpenalized.csv"
)
# Fall back to bordered txt if the wide csv is missing.
fcir_m_txt = file.path(
  project_root, "Simulation_Results", "FCIR_M", "tables",
  "New_FCIR_M_adaptive_lasso_selection_wide_unpenalized.txt"
)

out_dir = file.path(
  project_root, "Simulation_Results", "comparison_plots",
  "Ising_Theta_vs_FCIR_M_Theta_int"
)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Parse bordered ASCII tables written by the evaluate_* scripts.
read_bordered_table <- function(path) {
  lines = readLines(path, warn = FALSE)
  lines = lines[grepl("^\\|", lines)]
  # Drop separator-only rows that contain only | + - + spaces
  lines = lines[!grepl("^\\|[-+| ]+\\|$", lines)]
  parsed = lapply(lines, function(ln) {
    parts = strsplit(ln, "\\|", fixed = FALSE)[[1]]
    trimws(parts[nzchar(trimws(parts))])
  })
  header = parsed[[1]]
  body = do.call(rbind, parsed[-1])
  df = as.data.frame(body, stringsAsFactors = FALSE)
  names(df) = header
  num_cols = setdiff(header, "Matrix")
  df[num_cols] = lapply(df[num_cols], function(x) as.numeric(x))
  df
}

ising = read_bordered_table(ising_txt) %>%
  filter(Matrix == "Theta") %>%
  mutate(Model = "Ising (Theta)")

if (file.exists(fcir_m_csv)) {
  fcir_m = read.csv(fcir_m_csv, check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(Matrix == "Theta_int") %>%
    mutate(Model = "FCIR_M (Theta_int)")
} else {
  fcir_m = read_bordered_table(fcir_m_txt) %>%
    filter(Matrix == "Theta_int") %>%
    mutate(Model = "FCIR_M (Theta_int)")
}

# Harmonize column names (txt / csv may differ slightly).
rename_metrics <- function(df) {
  nm = names(df)
  nm[nm %in% c("no of correct zeros", "no.of.correct.zeros")] = "correct_zeros"
  nm[nm %in% c("no of correct non-zeros", "no.of.correct.non.zeros")] = "correct_nonzeros"
  names(df) = nm
  df
}

ising = rename_metrics(ising)
fcir_m = rename_metrics(fcir_m)

keep = c("Model", "N", "P", "n_nonzero", "n_total",
         "correct_zeros", "correct_nonzeros", "F1", "RMSE")
cmp = bind_rows(ising[, keep], fcir_m[, keep]) %>%
  mutate(
    P_lab = paste0("P = ", P),
    n_true_zeros = n_total - n_nonzero,
    zero_recovery = correct_zeros / n_true_zeros,
    nonzero_recovery = correct_nonzeros / n_nonzero
  )

# Sanity: only keep (N,P) present in both models.
common = cmp %>%
  count(N, P, name = "n_models") %>%
  filter(n_models == 2) %>%
  select(N, P)
cmp = semi_join(cmp, common, by = c("N", "P"))

write.csv(cmp, file.path(out_dir, "Theta_vs_Theta_int_metrics.csv"),
          row.names = FALSE)

metric_long = cmp %>%
  select(Model, N, P_lab,
         `Correct zeros` = correct_zeros,
         `Correct nonzeros` = correct_nonzeros,
         F1, RMSE,
         `Zero recovery rate` = zero_recovery,
         `Nonzero recovery rate` = nonzero_recovery) %>%
  pivot_longer(cols = -c(Model, N, P_lab),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(
    Metric,
    levels = c("Correct zeros", "Correct nonzeros", "F1", "RMSE",
               "Zero recovery rate", "Nonzero recovery rate")
  ))

# Main comparison: absolute counts + F1 + RMSE (what the tables report).
plot_main = metric_long %>%
  filter(Metric %in% c("Correct zeros", "Correct nonzeros", "F1", "RMSE")) %>%
  ggplot(aes(x = N, y = Value, colour = Model, shape = Model, group = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_continuous(trans = "log2", breaks = sort(unique(cmp$N))) +
  facet_grid(Metric ~ P_lab, scales = "free_y") +
  labs(
    title = "Ising Theta vs FCIR_M Theta_int (adaptive lasso)",
    x = "N (log2 scale)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "compare_Theta_vs_Theta_int_main.png"),
       plot_main, width = 9, height = 10, dpi = 150)

# Rate version: more fair when interpreting selection quality.
plot_rates = metric_long %>%
  filter(Metric %in% c("Zero recovery rate", "Nonzero recovery rate", "F1", "RMSE")) %>%
  ggplot(aes(x = N, y = Value, colour = Model, shape = Model, group = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_continuous(trans = "log2", breaks = sort(unique(cmp$N))) +
  facet_grid(Metric ~ P_lab, scales = "free_y") +
  labs(
    title = "Ising Theta vs FCIR_M Theta_int (rates + F1 + RMSE)",
    x = "N (log2 scale)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "compare_Theta_vs_Theta_int_rates.png"),
       plot_rates, width = 9, height = 10, dpi = 150)

message("Saved comparison plots to: ", out_dir)
message("Matched (N, P) cells: ", nrow(common))
print(cmp %>% arrange(P, N, Model) %>%
        select(Model, N, P, correct_zeros, correct_nonzeros, F1, RMSE))
