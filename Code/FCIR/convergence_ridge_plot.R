library(ggplot2)
library(tidyr)
library(dplyr)

project_root = "E:/RStudio/thesis"
rdata_dir = file.path(project_root, "Simulation_Results", "Rdata_Sparse")
output_dir = file.path(project_root, "Simulation_Results", "plots", "boxplots_perparameter_ridge_sparse")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)
ridge_lambdas = c(0.5, 1, 2)

format_lambda_label <- function(lambda_value) {
  gsub("\\.", "p", as.character(lambda_value))
}

safe_parameter_name <- function(param) {
  safe_name = gsub("\n.*", "", param)
  safe_name = gsub("\\[", "_", safe_name)
  safe_name = gsub("\\]", "", safe_name)
  safe_name = gsub(",", "_", safe_name)
  safe_name
}

df_list = list()

for (lambda_value in ridge_lambdas) {
  lambda_label = format_lambda_label(lambda_value)
  
  for (n in Ns) {
    for (p in Ps) {
      filename = file.path(
        rdata_dir,
        paste0("FCIR_estimates_ridge_lambda", lambda_label, "_N", n, "_P", p, ".Rdata")
      )
      
      if (!file.exists(filename)) {
        warning(paste("File not found:", filename, "- Skipping."))
        next
      }
      
      env = new.env()
      load(filename, envir = env)
      
      B_reps = nrow(env$est_beta_0)
      L = env$L
      K = env$K
      
      # beta_0 and alpha_0 use relative error.
      err_beta = as.data.frame(sweep(sweep(env$est_beta_0, 2, env$beta_0, "-"), 2, abs(env$beta_0), "/"))
      colnames(err_beta) = paste0("beta_0[", 1:L, "]\n(True: ", round(env$beta_0, 3), ")")
      df_beta = pivot_longer(err_beta, cols = everything(), names_to = "Parameter", values_to = "Error")
      
      err_alpha = as.data.frame(sweep(sweep(env$est_alpha_0, 2, env$alpha_0, "-"), 2, abs(env$alpha_0), "/"))
      colnames(err_alpha) = paste0("alpha_0[", 1:L, "]\n(True: ", round(env$alpha_0, 3), ")")
      df_alpha = pivot_longer(err_alpha, cols = everything(), names_to = "Parameter", values_to = "Error")
      
      # Sparse B_mat and A_mat may contain true zeros, so use absolute error.
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
      
      df_temp = bind_rows(df_beta, df_alpha, df_B, df_A)
      df_temp$N = n
      df_temp$P = p
      df_temp$Lambda = lambda_value
      df_temp$N_factor = factor(n, levels = Ns)
      df_temp$P_label = paste0("Species Size (P = ", p, ")")
      df_temp$Lambda_label = paste0("Ridge lambda = ", lambda_value)
      
      df_list[[paste(lambda_label, n, p, sep = "_")]] = df_temp
      message("Loaded ridge estimates: lambda=", lambda_value, ", N=", n, ", P=", p)
    }
  }
}

if (length(df_list) == 0) {
  stop("No ridge estimate files were loaded.")
}

df_all = bind_rows(df_list)
df_all$N_factor = factor(df_all$N, levels = Ns)
df_all$P_label = factor(df_all$P_label, levels = paste0("Species Size (P = ", Ps, ")"))
df_all$Lambda_label = factor(df_all$Lambda_label, levels = paste0("Ridge lambda = ", ridge_lambdas))

parameters = unique(df_all$Parameter)

for (param in parameters) {
  df_param = df_all %>% filter(.data$Parameter == param)
  
  y_label = if (grepl("^(beta|alpha)", param)) {
    "Relative Error: (Estimated - True) / |True|"
  } else {
    "Absolute Error: (Estimated - True)"
  }
  
  plot_obj = ggplot(df_param, aes(x = .data$N_factor, y = .data$Error)) +
    geom_boxplot(fill = "skyblue", color = "black",
                 outlier.shape = 16, outlier.alpha = 0.3, width = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    facet_grid(rows = vars(.data$P_label), cols = vars(.data$Lambda_label)) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Ridge Convergence of Parameter:", param),
         x = "Sample Size (N)",
         y = y_label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "lightgray", color = NA),
          strip.text = element_text(face = "bold", size = 12))
  
  plot_filename = file.path(
    output_dir,
    paste0("esterror_samplesize_ridge_", safe_parameter_name(param), ".png")
  )
  
  ggsave(filename = plot_filename, plot = plot_obj, width = 14, height = 8, dpi = 300)
  message("Saved plot: ", plot_filename)
}

message("Finished ridge per-parameter convergence plots.")
