library(ggplot2)
library(tidyr)
library(dplyr)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

for (n in Ns) {
  for (p in Ps) {
    
    filename = paste0("Simulation_Results/FCIR_estimates_N", n, "_P", p, ".Rdata")
    if (!file.exists(filename)) {
      warning(paste("File not found:", filename, "- Skipping."))
      next 
    }
    
    load(filename)
    print(paste("Successfully loaded combination: N =", n, ", P =", p))
    # 1. Beta_0
    err_beta = as.data.frame(sweep(est_beta_0, 2, beta_0, "-"))
    colnames(err_beta) = paste0("beta_0[", 1:ncol(est_beta_0), "]")
    df_beta = pivot_longer(err_beta, cols = everything(), names_to = "Parameter", values_to = "Error")
    mse_beta = mean(df_beta$Error^2, na.rm = TRUE)
    df_beta$Group = paste0("1. Env Main Effects (beta_0)\nOverall MSE: ", sprintf("%.4f", mse_beta))
    # 2. Alpha_0
    err_alpha = as.data.frame(sweep(est_alpha_0, 2, alpha_0, "-"))
    colnames(err_alpha) = paste0("alpha_0[", 1:ncol(est_alpha_0), "]")
    df_alpha = pivot_longer(err_alpha, cols = everything(), names_to = "Parameter", values_to = "Error")
    mse_alpha = mean(df_alpha$Error^2, na.rm = TRUE)
    df_alpha$Group = paste0("2. Trait Main Effects (alpha_0)\nOverall MSE: ", sprintf("%.4f", mse_alpha))
    # 3. B_mat (Environment-Trait Interactions)
    err_B = matrix(NA, nrow = B_reps, ncol = L * K)
    col_names_B = c()
    idx = 1
    for (l in 1:L) {
      for (k in 1:K) {
        err_B[, idx] = est_B_mat[l, k, ] - B_mat[l, k]
        col_names_B = c(col_names_B, paste0("B[", l, ",", k, "]"))
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
        col_names_A = c(col_names_A, paste0("A[", i, ",", j, "]"))
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
      facet_wrap(~ Group, scales = "free_x", ncol = 2) +
      theme_minimal(base_size = 14) +
      labs(title = paste("Estimation Error of ALL Parameters (N =", n, ", P =", p, ")"),
           subtitle = paste("Based on", B_reps, "Monte Carlo replications"),
           x = "Parameters",
           y = "Estimation Error (Estimated - True)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold", size = 12))
    
    plot_filename = paste0("Simulation_Results/Boxplot_ALL_Error_N", n, "_P", p, ".png")
    ggsave(filename = plot_filename, plot = plot_obj, width = 12, height = 9, dpi = 300)
    
    print(paste("--> Comprehensive Plot saved:", plot_filename))
  }
}

print("Finish！")

