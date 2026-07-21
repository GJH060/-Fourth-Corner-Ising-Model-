library(ggplot2)
library(tidyr)
library(dplyr)

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-/"
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Sparse")
plot_dir = file.path(project_root, "Simulation_Results", "FCIR", "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)
ridge_lambdas = c(5e-3, 5e-2, 5e-1)

format_lambda_label <- function(lambda_value) {
  gsub("\\.", "p", as.character(lambda_value))
}

ridge_methods = paste0("ridge_lambda", sapply(ridge_lambdas, format_lambda_label))
# Unpenalized all-parameter plots live in boxplots_allparam_unpenalized_sparse;
# this script only produces the fixed-lambda ridge all-parameter plots.
methods = "unpenalized"

for (method in methods) {
  for (n in Ns) {
    for (p in Ps) {
      
      if (method == "unpenalized") {
        filename = file.path(rdata_dir, paste0("FCIR_estimates_N", n, "_P", p, ".Rdata"))
      } else if (method == "penalized") {
        filename = file.path(rdata_dir, paste0("FCIR_estimates_penalized_N", n, "_P", p, ".Rdata"))
      } else if (grepl("^ridge_lambda", method)) {
        filename = file.path(rdata_dir, paste0("FCIR_estimates_", method, "_N", n, "_P", p, ".Rdata"))
      } else {
        warning(paste("Unknown method:", method, "- Skipping."))
        next
      }
    if (!file.exists(filename)) {
      warning(paste("File not found:", filename, "- Skipping."))
      next 
    }
    
    load(filename)
    print(paste("Successfully loaded combination: N =", n, ", P =", p))
    # 1. Beta_0
    err_beta = as.data.frame(sweep(sweep(est_beta_0, 2, beta_0, "-"), 2, abs(beta_0), "/"))
    colnames(err_beta) = paste0("beta_0[", seq_len(ncol(est_beta_0)), "]\n(True: ", round(beta_0, 3), ")")
    df_beta = pivot_longer(err_beta, cols = everything(), names_to = "Parameter", values_to = "Error")
    re_beta = mean(abs(df_beta$Error), na.rm = TRUE)
    df_beta$Group = paste0("1. Env Main Effects (beta_0)\nOverall MRE: ", sprintf("%.4f", re_beta))
    
    # 2. Alpha_0
    err_alpha = as.data.frame(sweep(sweep(est_alpha_0, 2, alpha_0, "-"), 2, abs(alpha_0), "/"))
    colnames(err_alpha) = paste0("alpha_0[", seq_len(ncol(est_alpha_0)), "]\n(True: ", round(alpha_0, 3), ")")
    df_alpha = pivot_longer(err_alpha, cols = everything(), names_to = "Parameter", values_to = "Error")
    re_alpha = mean(abs(df_alpha$Error), na.rm = TRUE)
    df_alpha$Group = paste0("2. Trait Main Effects (alpha_0)\nOverall MRE: ", sprintf("%.4f", re_alpha))
    
    # 3. B_mat (Environment-Trait Interactions)
    err_B = matrix(NA, nrow = B_reps, ncol = L * K)
    col_names_B = c()
    idx = 1
    for (l in 1:L) {
      for (k in 1:K) {
        err_B[, idx] = est_B_mat[l, k, ] - B_mat[l, k]
        col_names_B = c(col_names_B, paste0("B[", l, ",", k, "]\n(True: ", round(B_mat[l, k], 3), ")"))
        idx = idx + 1
      }
    }
    err_B = as.data.frame(err_B)
    colnames(err_B) = col_names_B
    df_B = pivot_longer(err_B, cols = everything(), names_to = "Parameter", values_to = "Error")
    
    mse_B = mean(df_B$Error^2, na.rm = TRUE)
    df_B$Group = paste0("3. Env-Trait Interactions (B_mat)\nOverall MSE: ", sprintf("%.4f", mse_B))
    
    # 4. Matrix A (A_mat)
    dim_A1 = dim(est_A_mat)[1]
    dim_A2 = dim(est_A_mat)[2]
    err_A = matrix(NA, nrow = B_reps, ncol = dim_A1 * dim_A2)
    col_names_A = c()
    idx = 1
    for (i in 1:dim_A1) {
      for (j in 1:dim_A2) {
        err_A[, idx] = est_A_mat[i, j, ] - A_mat[i, j]
        col_names_A = c(col_names_A, paste0("A[", i, ",", j, "]\n(True: ", round(A_mat[i, j], 3), ")"))
        idx = idx + 1
      }
    }
    err_A = as.data.frame(err_A)
    colnames(err_A) = col_names_A
    df_A = pivot_longer(err_A, cols = everything(), names_to = "Parameter", values_to = "Error")
    mse_A = mean(df_A$Error^2, na.rm = TRUE)
    df_A$Group = paste0("4. Matrix A (A_mat)\nOverall MSE: ", sprintf("%.4f", mse_A))
    
    df_all = bind_rows(df_beta, df_alpha, df_B, df_A)
    
    plot_obj <- ggplot(df_all, aes(x = Parameter, y = Error)) +
      geom_boxplot(fill = "skyblue", color = "black", 
                   outlier.shape = 16, outlier.alpha = 0.4, width = 0.5) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~ Group, scales = "free", ncol = 2) +
      theme_minimal(base_size = 14) +
      labs(title = paste("Estimation Error (", method, ", N =", n, ", P =", p, ")"),
           subtitle = paste("Based on", B_reps, "Monte Carlo replications"),
           x = "Parameters",
           y = "Error (Rel. Error for beta/alpha, Abs. Error for B/A)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold", size = 12))
    
    mixederr_dir = file.path(plot_dir, "boxplots_allparam_ridge_fixedlambda_sparse")
    if (!dir.exists(mixederr_dir)) dir.create(mixederr_dir, recursive = TRUE)
    plot_filename = file.path(mixederr_dir, paste0("Boxplot_ALL_MixedError_", method, "_N", n, "_P", p, ".png"))
    ggsave(filename = plot_filename, plot = plot_obj, width = 12, height = 9, dpi = 300)
    
    print(paste("--> Comprehensive Plot saved:", plot_filename))
  }
}
}

print("Finished")
