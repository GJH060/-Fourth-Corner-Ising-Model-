# Compare two ridge cross-validation foldings on the same simulated data:
#   - "Random folds"       : FCIR_estimates_ridge_cv_N*_P*.Rdata          (baseline)
#   - "Site-grouped folds" : FCIR_estimates_ridge_cv_grouped_N*_P*.Rdata  (new)
# Both use the same datasets and model; only how cv.glmnet builds folds differs.

library(ggplot2)
library(tidyr)
library(dplyr)

project_root = "E:/RStudio/thesis"
rdata_dir = file.path(project_root, "Simulation_Results", "Rdata_Sparse")
output_dir = file.path(project_root, "Simulation_Results", "plots", "boxplots_perparameter_ridge_cv_compare_sparse")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

# File infix -> readable label (order sets the legend / dodge order)
cv_variants = c(
  "ridge_cv"         = "Random folds",
  "ridge_cv_grouped" = "Site-grouped folds"
)
fill_values = c("Random folds" = "skyblue", "Site-grouped folds" = "lightgreen")

safe_parameter_name <- function(param) {
  safe_name = gsub("\n.*", "", param)
  safe_name = gsub("\\[", "_", safe_name)
  safe_name = gsub("\\]", "", safe_name)
  safe_name = gsub(",", "_", safe_name)
  safe_name
}

# Per-parameter long error data frame for one loaded estimate file.
# beta_0/alpha_0 use relative error; sparse B/A use absolute error.
build_error_long <- function(env) {
  B_reps = nrow(env$est_beta_0)
  L = env$L
  K = env$K

  err_beta = as.data.frame(sweep(sweep(env$est_beta_0, 2, env$beta_0, "-"), 2, abs(env$beta_0), "/"))
  colnames(err_beta) = paste0("beta_0[", 1:L, "]\n(True: ", round(env$beta_0, 3), ")")
  df_beta = pivot_longer(err_beta, cols = everything(), names_to = "Parameter", values_to = "Error")

  err_alpha = as.data.frame(sweep(sweep(env$est_alpha_0, 2, env$alpha_0, "-"), 2, abs(env$alpha_0), "/"))
  colnames(err_alpha) = paste0("alpha_0[", 1:L, "]\n(True: ", round(env$alpha_0, 3), ")")
  df_alpha = pivot_longer(err_alpha, cols = everything(), names_to = "Parameter", values_to = "Error")

  err_B = matrix(NA, nrow = B_reps, ncol = L * K)
  col_names_B = c()
  idx = 1
  for (l in 1:L) {
    for (k in 1:K) {
      err_B[, idx] = env$est_B_mat[l, k, ] - env$B_mat[l, k]
      col_names_B = c(col_names_B, paste0("B[", l, ",", k, "]\n(True: ", round(env$B_mat[l, k], 3), ")"))
      idx = idx + 1
    }
  }
  err_B = as.data.frame(err_B)
  colnames(err_B) = col_names_B
  df_B = pivot_longer(err_B, cols = everything(), names_to = "Parameter", values_to = "Error")

  err_A = matrix(NA, nrow = B_reps, ncol = L * K)
  col_names_A = c()
  idx = 1
  for (l in 1:L) {
    for (k in 1:K) {
      err_A[, idx] = env$est_A_mat[l, k, ] - env$A_mat[l, k]
      col_names_A = c(col_names_A, paste0("A[", l, ",", k, "]\n(True: ", round(env$A_mat[l, k], 3), ")"))
      idx = idx + 1
    }
  }
  err_A = as.data.frame(err_A)
  colnames(err_A) = col_names_A
  df_A = pivot_longer(err_A, cols = everything(), names_to = "Parameter", values_to = "Error")

  bind_rows(df_beta, df_alpha, df_B, df_A)
}

df_list = list()
lambda_list = list()

for (infix in names(cv_variants)) {
  variant_label = cv_variants[[infix]]

  for (n in Ns) {
    for (p in Ps) {
      filename = file.path(
        rdata_dir,
        paste0("FCIR_estimates_", infix, "_N", n, "_P", p, ".Rdata")
      )

      if (!file.exists(filename)) {
        warning(paste("File not found:", filename, "- Skipping."))
        next
      }

      env = new.env()
      load(filename, envir = env)

      df = build_error_long(env)
      df$N = n
      df$P = p
      df$CV_type = variant_label
      df_list[[paste(infix, n, p, sep = "_")]] = df

      if (!is.null(env$est_lambda)) {
        lambda_list[[paste(infix, n, p, sep = "_")]] = data.frame(
          N = n, P = p, CV_type = variant_label, lambda = env$est_lambda
        )
      }

      message("Loaded ", variant_label, ": N=", n, ", P=", p)
    }
  }
}

if (length(df_list) == 0) {
  stop("No ridge CV estimate files were loaded.")
}

cv_levels = unname(cv_variants)

df_all = bind_rows(df_list)
df_all$N_factor = factor(df_all$N, levels = Ns)
df_all$P_label = factor(paste0("Species Size (P = ", df_all$P, ")"),
                        levels = paste0("Species Size (P = ", Ps, ")"))
df_all$CV_type = factor(df_all$CV_type, levels = cv_levels)

loaded_levels = intersect(cv_levels, unique(as.character(df_all$CV_type)))
if (length(loaded_levels) < 2) {
  warning(paste0("Only one CV variant was found (", paste(loaded_levels, collapse = ", "),
                 "). Run main.r to generate the site-grouped results before comparing."))
}

parameters = unique(df_all$Parameter)

for (param in parameters) {
  df_param = df_all %>% filter(.data$Parameter == param)

  y_label = if (grepl("^(beta|alpha)", param)) {
    "Relative Error: (Estimated - True) / |True|"
  } else {
    "Absolute Error: (Estimated - True)"
  }

  plot_obj = ggplot(df_param, aes(x = .data$N_factor, y = .data$Error, fill = .data$CV_type)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    geom_boxplot(color = "black", outlier.shape = 16, outlier.alpha = 0.3,
                 position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ .data$P_label, ncol = 2) +
    scale_fill_manual(values = fill_values) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Ridge CV folding comparison of Parameter:", param),
         x = "Sample Size (N)",
         y = y_label,
         fill = "CV fold type") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top",
          strip.background = element_rect(fill = "lightgray", color = NA),
          strip.text = element_text(face = "bold", size = 12))

  plot_filename = file.path(
    output_dir,
    paste0("compare_esterror_samplesize_", safe_parameter_name(param), ".png")
  )

  ggsave(filename = plot_filename, plot = plot_obj, width = 13, height = 6, dpi = 300)
  message("Saved plot: ", plot_filename)
}

# --- CV-selected lambda: random vs site-grouped folds ---
if (length(lambda_list) > 0) {
  df_lambda = bind_rows(lambda_list)
  df_lambda$N_factor = factor(df_lambda$N, levels = Ns)
  df_lambda$P_label = factor(paste0("Species Size (P = ", df_lambda$P, ")"),
                             levels = paste0("Species Size (P = ", Ps, ")"))
  df_lambda$CV_type = factor(df_lambda$CV_type, levels = cv_levels)

  lambda_plot = ggplot(df_lambda, aes(x = .data$N_factor, y = .data$lambda, fill = .data$CV_type)) +
    geom_boxplot(color = "black", outlier.shape = 16, outlier.alpha = 0.3,
                 position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ .data$P_label, ncol = 2) +
    scale_fill_manual(values = fill_values) +
    theme_minimal(base_size = 14) +
    labs(title = "CV-selected Ridge lambda.min: random vs site-grouped folds",
         x = "Sample Size (N)",
         y = "Selected lambda",
         fill = "CV fold type") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top",
          strip.background = element_rect(fill = "lightgray", color = NA),
          strip.text = element_text(face = "bold", size = 12))

  lambda_filename = file.path(output_dir, "compare_selected_lambda.png")
  ggsave(filename = lambda_filename, plot = lambda_plot, width = 13, height = 6, dpi = 300)
  message("Saved plot: ", lambda_filename)
}

message("Finished ridge CV folding comparison plots.")
