library(ggplot2)
library(tidyr)
library(dplyr)

generate_parameter_convergence_plots <- function(Ns = c(50, 100, 200, 400, 800),
                                                 Ps = c(30, 60),
                                                 L = 3, K = 2,
                                                 B_reps = 1000,
                                                 result_dir = "Simulation_Results",
                                                 output_dir = "Convergence_Plots") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  df_list <- list() 
  
  for (n in Ns) {
    for (p in Ps) {
      
      filename = paste0(result_dir, "/FCIR_estimates_N", n, "_P", p, ".Rdata")
      if (!file.exists(filename)) {
        warning(paste())
        next 
      }
      
      load(filename)
      

      err_beta = as.data.frame(sweep(est_beta_0, 2, beta_0, "-"))
      colnames(err_beta) = paste0("beta_0[", 1:L, "]")
      df_beta = pivot_longer(err_beta, cols = everything(), names_to = "Parameter", values_to = "Error")
      

      err_alpha = as.data.frame(sweep(est_alpha_0, 2, alpha_0, "-"))
      colnames(err_alpha) = paste0("alpha_0[", 1:L, "]")
      df_alpha = pivot_longer(err_alpha, cols = everything(), names_to = "Parameter", values_to = "Error")
      

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
      

      err_A = matrix(NA, nrow = B_reps, ncol = L * K)
      col_names_A = c()
      idx = 1
      for (i in 1:L) {
        for (j in 1:K) {
          err_A[, idx] = est_A_mat[i, j, ] - A_mat[i, j]
          col_names_A = c(col_names_A, paste0("A[", i, ",", j, "]"))
          idx = idx + 1
        }
      }
      err_A = as.data.frame(err_A)
      colnames(err_A) = col_names_A
      df_A = pivot_longer(err_A, cols = everything(), names_to = "Parameter", values_to = "Error")
      

      df_temp = bind_rows(df_beta, df_alpha, df_B, df_A)
      df_temp$N = n
      df_temp$P = p
      
      df_list[[paste(n, p)]] = df_temp
    }
  }
  
 
  df_all = bind_rows(df_list)
  
  df_all$N_factor = factor(df_all$N, levels = Ns)

  df_all$P_label = paste0("Species Size (P = ", df_all$P, ")")
  

  parameters = unique(df_all$Parameter)
  
  
  for (param in parameters) {
    
    df_param = df_all %>% filter(Parameter == param)
    
    plot_obj <- ggplot(df_param, aes(x = N_factor, y = Error)) +
      geom_boxplot(fill = "skyblue", color = "black", 
                   outlier.shape = 16, outlier.alpha = 0.3, width = 0.5) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
      facet_wrap(~ P_label, ncol = 2) +
      theme_minimal(base_size = 14) +
      labs(title = paste("Convergence of Parameter:", param),
           x = "Sample Size (N)",
           y = "Estimation Error (Estimated - True)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold", size = 12))
    
    safe_param_name = gsub("\\[", "_", param)
    safe_param_name = gsub("\\]", "", safe_param_name)
    safe_param_name = gsub(",", "_", safe_param_name)
    
    plot_filename = paste0(output_dir, "/Convergence_", safe_param_name, ".png")
    
    ggsave(filename = plot_filename, plot = plot_obj, width = 10, height = 6, dpi = 300)
    
  }
  
}

generate_parameter_convergence_plots()