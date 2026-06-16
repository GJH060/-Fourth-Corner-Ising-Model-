# Mean signed interaction-edge error for Theta_s[j, j'] (j < j')
#
# In generate_FCIR.r, for site s and species pair (j, j'):
#   theta_{jj'}(s) = sum(x_s * alpha_0) + sum(x_s * (A_mat %*% Delta[j, j', ]))
# Pseudo-likelihood gives (hat_alpha_0, hat_A_mat); plug in the same Delta and X
# to obtain hat_theta_{jj'}(s).
#
# We report the average SIGNED error over all sites and unordered pairs:
#   (2 / (N * P * (P - 1))) * sum_s sum_{j<j'} (hat_theta - theta_true)
# (This is the natural "one number per replicate" analogue of an edge-mean;
#  1/(N*P) would not match the number of edges unless redefined.)

library(ggplot2)
library(dplyr)

# --- paths (edit if your project root differs) ---
project_root <- "E:/RStudio/thesis"
# Sparse (matches generate_FCIR.r / main.r outputs)
data_dir <- file.path(project_root, "Simulation_Results", "Rdata_Sparse")
est_dir <- file.path(project_root, "Simulation_Results", "Rdata_Sparse")
# Dense workflow: set use_dense <- TRUE and point data_dir / est_dir to Rdata_Dense.
use_dense <- FALSE
if (use_dense) {
  data_dir <- file.path(project_root, "Simulation_Results", "Rdata_Dense")
  est_dir <- file.path(project_root, "Simulation_Results", "Rdata_Dense")
}
scenario_tag <- if (use_dense) "Dense" else "Sparse"

plot_dir <- file.path(project_root, "Simulation_Results", "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

Ns <- c(50, 100, 200, 400, 800)
Ps <- c(30, 60)

ridge_lambdas <- c(0.5, 1, 2)

format_lambda_label <- function(lambda_value) {
  gsub("\\.", "p", as.character(lambda_value))
}

ridge_methods <- paste0("ridge_lambda", sapply(ridge_lambdas, format_lambda_label))
methods <- c("unpenalized", ridge_methods)

#' One column per unordered pair (j, j'), j < j': trait difference vector in R^K
pair_trait_diff_matrix <- function(Tr) {
  P <- nrow(Tr)
  K <- ncol(Tr)
  prs <- combn(P, 2)
  M <- ncol(prs)
  Dmat <- matrix(NA_real_, nrow = K, ncol = M)
  for (m in seq_len(M)) {
    j <- prs[1, m]
    jp <- prs[2, m]
    Dmat[, m] <- abs(Tr[j, ] - Tr[jp, ])
  }
  list(Dmat = Dmat, M = M, prs = prs)
}

#' For one Monte Carlo replicate b: scalar mean signed error over edges x sites
mean_theta_edge_signed_error <- function(X, Tr, alpha_0, A_mat, est_alpha_row, est_A_mat) {
  N <- nrow(X)
  pr <- pair_trait_diff_matrix(Tr)
  Dmat <- pr$Dmat
  M <- pr$M

  # True contributions from interactions: for each site, length-M vector
  err_sum <- 0.0
  n_terms <- N * M

  A_true <- A_mat
  a0_true <- alpha_0
  A_hat <- est_A_mat
  a0_hat <- est_alpha_row

  W_true <- A_true %*% Dmat # L x M
  W_hat <- A_hat %*% Dmat

  for (s in seq_len(N)) {
    x_s <- X[s, ]
    base_true <- sum(x_s * a0_true)
    base_hat <- sum(x_s * a0_hat)
    sec_true <- colSums(sweep(W_true, 1L, x_s, `*`))
    sec_hat <- colSums(sweep(W_hat, 1L, x_s, `*`))
    theta_true <- base_true + sec_true
    theta_hat <- base_hat + sec_hat
    err_sum <- err_sum + sum(theta_hat - theta_true)
  }

  err_sum / n_terms
}

for (method in methods) {
  out1 <- file.path(
    plot_dir,
    paste0("Theta_edge_mean_signed_error_boxplot_", scenario_tag, "_", method, ".png")
  )
  out2 <- file.path(
    plot_dir,
    paste0("Theta_edge_mean_signed_error_line_", scenario_tag, "_", method, ".png")
  )
  
  if (file.exists(out1) && file.exists(out2)) {
    message("Plots already exist for method=", method, ", scenario=", scenario_tag, " - Skipping.")
    next
  }
  
  results <- list()
  
  for (n in Ns) {
    for (p in Ps) {
      if (use_dense) {
        data_file <- file.path(data_dir, sprintf("FCIR_data_Dense_N%d_P%d.Rdata", n, p))
        if (method == "unpenalized") {
          est_file <- file.path(est_dir, sprintf("FCIR_estimates_Dense_N%d_P%d.Rdata", n, p))
        } else if (grepl("^ridge_lambda", method)) {
          est_file <- file.path(est_dir, sprintf("FCIR_estimates_%s_Dense_N%d_P%d.Rdata", method, n, p))
        } else {
          est_file <- file.path(est_dir, sprintf("FCIR_estimates_penalized_Dense_N%d_P%d.Rdata", n, p))
        }
      } else {
        data_file <- file.path(data_dir, sprintf("FCIR_data_N%d_P%d.Rdata", n, p))
        if (method == "unpenalized") {
          est_file <- file.path(est_dir, sprintf("FCIR_estimates_N%d_P%d.Rdata", n, p))
        } else if (grepl("^ridge_lambda", method)) {
          est_file <- file.path(est_dir, sprintf("FCIR_estimates_%s_N%d_P%d.Rdata", method, n, p))
        } else {
          est_file <- file.path(est_dir, sprintf("FCIR_estimates_penalized_N%d_P%d.Rdata", n, p))
        }
      }
  
      if (!file.exists(data_file) || !file.exists(est_file)) {
        message("Skip N=", n, " P=", p, " method=", method, " (missing data or estimates)")
        next
      }
  
      load(data_file)
      load(est_file)
  
      B_reps <- dim(est_A_mat)[3]
      vals <- numeric(B_reps)
      for (b in seq_len(B_reps)) {
        vals[b] <- mean_theta_edge_signed_error(
          X, Tr, alpha_0, A_mat,
          est_alpha_0[b, ], est_A_mat[, , b]
        )
      }
  
      results[[paste0("N", n, "_P", p)]] <- data.frame(
        N = n, P = p, B = seq_len(B_reps),
        mean_signed_edge_error = vals
      )
      message("Done N=", n, " P=", p, " method=", method)
    }
  }
  
  if (length(results) == 0) {
    warning(paste("No estimate files loaded for method:", method))
    next
  }
  
  df_all <- bind_rows(results)
  df_all$N_f <- factor(df_all$N, levels = Ns)
  df_all$P_f <- factor(df_all$P, levels = Ps)
  
  # --- Plot 1: distribution across Monte Carlo replicates ---
  p1 <- ggplot(df_all, aes(x = N_f, y = mean_signed_edge_error, fill = P_f)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(outlier.alpha = 0.35, position = position_dodge(width = 0.85)) +
    labs(
      title = paste0(
        "Mean signed error of Ising pair potentials Theta_s[j,j'] (",
        method, ", ", scenario_tag, ")"
      ),
      subtitle = paste0(
        "For each replicate: average of (hat(theta) - theta) over s=1..N and all unordered pairs j<j'."
      ),
      x = "Sample size N",
      y = "Mean signed edge error",
      fill = "P"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  
  ggsave(out1, p1, width = 9, height = 6, dpi = 150)
  message("Saved: ", out1)
  
  # --- Plot 2: Monte Carlo mean +/- sd vs N (line) ---
  summ <- df_all %>%
    group_by(N, P) %>%
    summarise(
      mc_mean = mean(mean_signed_edge_error, na.rm = TRUE),
      mc_sd = stats::sd(mean_signed_edge_error, na.rm = TRUE),
      .groups = "drop"
    )
  
  p2 <- ggplot(summ, aes(x = N, y = mc_mean, color = factor(P), group = factor(P))) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(ymin = mc_mean - mc_sd, ymax = mc_mean + mc_sd),
      width = 30,
      linewidth = 0.35
    ) +
    labs(
      title = paste0(
        "Monte Carlo summary of mean signed Theta edge error (",
        method, ", ", scenario_tag, ")"
      ),
      x = "Sample size N",
      y = "Mean across replicates of signed edge error",
      color = "P"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  
  ggsave(out2, p2, width = 9, height = 6, dpi = 150)
  message("Saved: ", out2)
}

message("Finished theta edge mean signed-error plots.")
