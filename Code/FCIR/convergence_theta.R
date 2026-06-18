
library(ggplot2)
library(dplyr)

project_root <- "E:/RStudio/thesis"
plot_dir <- file.path(project_root, "Simulation_Results", "plots", "theta_convergence")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

Ns <- c(50, 100, 200, 400, 800)
Ps <- c(30, 60)

# ---- methods: label -> file infix builder ----
method_levels <- c("Unpenalized", "Ridge lambda=0.005", "Ridge lambda=0.05",
                   "Ridge lambda=0.5", "Grouped-CV")
method_fill <- c("Unpenalized" = "#9E9E9E", "Ridge lambda=0.005" = "#A6CEE3",
                 "Ridge lambda=0.05" = "#5599C8", "Ridge lambda=0.5" = "#1F5C8B",
                 "Grouped-CV" = "#33A02C")

# estimate-file stem (without dir / "_N{n}_P{p}.Rdata"); infix = "" or "_Dense"
est_stem <- function(method, infix) {
  switch(method,
    "Unpenalized"      = paste0("FCIR_estimates", infix),
    "Ridge lambda=0.005" = paste0("FCIR_estimates_ridge_lambda0p005", infix),
    "Ridge lambda=0.05"  = paste0("FCIR_estimates_ridge_lambda0p05", infix),
    "Ridge lambda=0.5"   = paste0("FCIR_estimates_ridge_lambda0p5", infix),
    "Grouped-CV"       = paste0("FCIR_estimates_ridge_cv_grouped", infix)
  )
}

# K x M matrix of pairwise trait differences over unordered pairs j<j'
pair_delta <- function(Tr) {
  prs <- combn(nrow(Tr), 2)
  Dmat <- apply(prs, 2, function(pr) abs(Tr[pr[1], ] - Tr[pr[2], ]))
  matrix(Dmat, nrow = ncol(Tr))   # K x M
}

# per-replicate RMSE of theta from parameter differences
rmse_main <- function(X, TrT, dbeta0, dB) {
  # E[s,j] = x_s . dbeta0 + (X %*% dB %*% TrT)[s,j]
  E <- matrix(as.numeric(X %*% dbeta0), nrow = nrow(X), ncol = ncol(TrT)) + X %*% dB %*% TrT
  sqrt(mean(E^2))
}
rmse_edge <- function(X, Dmat, dalpha0, dA) {
  E <- matrix(as.numeric(X %*% dalpha0), nrow = nrow(X), ncol = ncol(Dmat)) + X %*% (dA %*% Dmat)
  sqrt(mean(E^2))
}

scenarios <- c("Dense", "Sparse")

for (scenario in scenarios) {
  infix   <- if (scenario == "Dense") "_Dense" else ""
  data_dir <- file.path(project_root, "Simulation_Results",
                        if (scenario == "Dense") "Rdata_Dense" else "Rdata_Sparse")

  rows_main <- list()
  rows_edge <- list()

  for (n in Ns) {
    for (p in Ps) {
      data_file <- file.path(data_dir, sprintf("FCIR_data%s_N%d_P%d.Rdata", infix, n, p))
      if (!file.exists(data_file)) { message("Missing data: ", data_file); next }
      d <- new.env(); load(data_file, envir = d)
      X <- d$X; Tr <- d$Tr
      TrT  <- t(Tr)              # K x P
      Dmat <- pair_delta(Tr)     # K x M

      for (method in method_levels) {
        est_file <- file.path(data_dir, sprintf("%s_N%d_P%d.Rdata", est_stem(method, infix), n, p))
        if (!file.exists(est_file)) { message("Missing est: ", basename(est_file)); next }
        e <- new.env(); load(est_file, envir = e)

        B_reps <- dim(e$est_A_mat)[3]
        vm <- numeric(B_reps); ve <- numeric(B_reps)
        for (b in seq_len(B_reps)) {
          vm[b] <- rmse_main(X, TrT,
                             e$est_beta_0[b, ] - e$beta_0,
                             e$est_B_mat[, , b] - e$B_mat)
          ve[b] <- rmse_edge(X, Dmat,
                             e$est_alpha_0[b, ] - e$alpha_0,
                             e$est_A_mat[, , b] - e$A_mat)
        }
        key <- sprintf("%s_N%d_P%d", method, n, p)
        rows_main[[key]] <- data.frame(N = n, P = p, method = method, rmse = vm)
        rows_edge[[key]] <- data.frame(N = n, P = p, method = method, rmse = ve)
        message(scenario, " | ", method, " N=", n, " P=", p)
      }
    }
  }

  make_plot <- function(rows, theta_label, file_tag) {
    if (length(rows) == 0) { warning("No data for ", scenario, " ", file_tag); return(invisible()) }
    df <- bind_rows(rows)
    df$N_f <- factor(df$N, levels = Ns)
    df$P_label <- factor(paste0("Species Size (P = ", df$P, ")"),
                         levels = paste0("Species Size (P = ", Ps, ")"))
    df$method <- factor(df$method, levels = method_levels)

    g <- ggplot(df, aes(x = N_f, y = rmse, fill = method)) +
      geom_boxplot(color = "black", outlier.shape = 16, outlier.alpha = 0.25,
                   position = position_dodge(width = 0.85), width = 0.75, linewidth = 0.3) +
      facet_wrap(~ P_label, ncol = 2, scales = "free_y") +
      scale_fill_manual(values = method_fill) +
      theme_minimal(base_size = 14) +
      labs(title = paste0("Convergence of ", theta_label, " [", scenario, "]"),
           subtitle = "Per-replicate RMSE over all sites x (species or unordered pairs)",
           x = "Sample Size (N)", y = expression(RMSE~of~theta), fill = "Method") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "top",
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold", size = 12))

    out <- file.path(plot_dir, sprintf("theta_%s_convergence_%s.png", file_tag, scenario))
    ggsave(out, g, width = 14, height = 6, dpi = 300)
    message("Saved: ", out)
  }

  make_plot(rows_main, "main-effect theta[s,jj]", "main")
  make_plot(rows_edge, "interaction theta[s,jj']", "edge")
}

message("Finished theta convergence plots.")
