# Figures comparing the adaptive lasso with the ordinary lasso (full FCIR) and
# the unpenalized MPLE (all three models). Reads the tables written by
# compare_estimators_RMSE.r.

library(ggplot2)
library(patchwork)

project_root = "D:/ising model thesis/-Fourth-Corner-Ising-Model-"
sim_dir = file.path(project_root, "Simulation_Results")

Ns = c(50, 100, 200, 400, 800)
est_labels = c(adaptive_lasso = "Adaptive lasso", lasso = "Lasso", unpenalized = "Unpenalized MPLE")
est_colours = c("Adaptive lasso" = "#0072B2", "Lasso" = "#009E73", "Unpenalized MPLE" = "#D55E00")
rmse_cap = 10   # mean RMSEs above this are clipped to the top edge (diverged fits)

base_theme = theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_blank(), axis.line = element_line(colour = "black"))

x_scale = scale_x_continuous(trans = "log2", breaks = Ns, name = "Number of sites, N")

prep = function(df) {
  df$Estimator = factor(est_labels[df$Estimator], levels = est_labels)
  df$Pf = factor(paste0("P = ", df$P), levels = paste0("P = ", sort(unique(df$P))))
  df$clipped = df$Mean_RMSE > rmse_cap
  df$y = pmin(df$Mean_RMSE, rmse_cap)
  df
}

rmse_panel = function(df, title) {
  ggplot(df, aes(N, y, colour = Estimator, linetype = Pf, group = interaction(Estimator, Pf))) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(shape = clipped), size = 2.4) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_colour_manual(values = est_colours, drop = TRUE) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_y_log10(name = "Mean RMSE") +
    x_scale + ggtitle(title) + base_theme
}

save_rmse_figure = function(model_dir, blocks, file, ncol) {
  df = prep(read.csv(file.path(sim_dir, model_dir, "tables", "estimator_RMSE_comparison.csv")))
  panels = lapply(names(blocks), function(b) rmse_panel(df[df$Block == b, ], blocks[[b]]))
  p = wrap_plots(panels, ncol = ncol) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  out = file.path(sim_dir, model_dir, "plots", file)
  ggsave(out, p, width = 4.2 * ncol, height = 3.6 * ceiling(length(blocks) / ncol) + 0.8, dpi = 300, bg = "white")
  message("Saved: ", out)
}

save_rmse_figure("FCIR", list(beta_0 = expression(bold(beta)[0]), B = expression(bold(B)),
                              alpha_0 = expression(bold(alpha)[0]), A = expression(bold(A))),
                 "FCIR_RMSE_estimator_comparison.png", ncol = 2)
save_rmse_figure("FCIR_I", list(Beta = "Species-specific main effects",
                                A = expression(bold(A)), alpha_0 = expression(bold(alpha)[0])),
                 "FCIR_I_RMSE_estimator_comparison.png", ncol = 3)
save_rmse_figure("FCIR_M", list(beta_0 = expression(bold(beta)[0]), B = expression(bold(B)),
                                Theta = expression(bold(Theta))),
                 "FCIR_M_RMSE_estimator_comparison.png", ncol = 3)

# ---- Ordinary vs adaptive lasso: support recovery for B and A (full FCIR) ----
ad = read.csv(file.path(sim_dir, "FCIR", "tables", "adaptive_lasso_selection_wide_unpenalized_Sparse.csv"),
              check.names = FALSE)
ad = data.frame(Estimator = "adaptive_lasso", Block = ad$Matrix, N = ad$N, P = ad$P,
                correct_zeros = ad$`no of correct zeros`, correct_nonzeros = ad$`no of correct non-zeros`,
                F1 = ad$F1)
la = read.csv(file.path(sim_dir, "FCIR", "tables", "lasso_selection_wide.csv"))
la = data.frame(Estimator = "lasso", la[, c("Block", "N", "P", "correct_zeros", "correct_nonzeros", "F1")])
sel = rbind(ad, la)
sel$Estimator = factor(est_labels[sel$Estimator], levels = est_labels[1:2])
sel$Pf = factor(paste0("P = ", sel$P))

metrics = list(correct_zeros = "Correct zeros", correct_nonzeros = "Correct nonzeros", F1 = quote(F[1]))
panels = list()
for (b in c("B", "A")) {
  for (m in names(metrics)) {
    d = sel[sel$Block == b, ]
    d$value = d[[m]]
    ymax = if (m == "F1") 1 else 3
    panels[[length(panels) + 1]] = ggplot(d, aes(N, value, colour = Estimator, linetype = Pf,
                                                  group = interaction(Estimator, Pf))) +
      geom_hline(yintercept = ymax, linetype = "dotted", colour = "grey40") +
      geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
      scale_colour_manual(values = est_colours[1:2]) +
      scale_linetype_manual(values = c("solid", "dashed")) +
      scale_y_continuous(limits = c(0, ymax), name = as.expression(metrics[[m]])) +
      x_scale + ggtitle(as.expression(bquote(bold(.(as.name(b))) * ": " * .(metrics[[m]])))) + base_theme
  }
}
p = wrap_plots(panels, ncol = 3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
out = file.path(sim_dir, "FCIR", "plots", "FCIR_lasso_vs_adaptive_selection.png")
ggsave(out, p, width = 12.6, height = 8, dpi = 300, bg = "white")
message("Saved: ", out)
