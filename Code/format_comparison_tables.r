# Writes LaTeX table bodies for the estimator comparisons in the thesis appendix
# from the CSVs produced by compare_estimators_RMSE.r.

project_root = "D:/ising model thesis/-Fourth-Corner-Ising-Model-"
sim_dir = file.path(project_root, "Simulation_Results")
out_file = file.path(sim_dir, "comparison_tables.tex")

fmt = function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("--")
    if (abs(v) >= 1000) {
      k = floor(log10(abs(v)))
      return(sprintf("$%.1f\\times10^{%d}$", v / 10^k, k))
    }
    sprintf("$%.3f$", v)
  }, character(1))
}

rmse_rows = function(df, blocks, columns, p_first = FALSE) {
  lines = character(0)
  for (b in names(blocks)) {
    d = df[df$Block == b, ]
    keys = unique(d[if (p_first) order(d$P, d$N) else order(d$N, d$P), c("N", "P")])
    for (i in seq_len(nrow(keys))) {
      r = d[d$N == keys$N[i] & d$P == keys$P[i], ]
      vals = vapply(columns, function(cl) {
        x = r[r$Estimator == cl$est, cl$stat]
        if (length(x) == 0) NA_real_ else x
      }, numeric(1))
      lead = if (i == 1) sprintf("\\multirow{%d}{*}{%s}", nrow(keys), blocks[[b]]) else ""
      lines = c(lines, sprintf("%s & %d & %d & %s \\\\", lead, keys$N[i], keys$P[i],
                               paste(fmt(vals), collapse = " & ")))
    }
    lines = c(lines, "\\midrule")
  }
  head(lines, -1)
}

out = character(0)

fcir = read.csv(file.path(sim_dir, "FCIR", "tables", "estimator_RMSE_comparison.csv"))
out = c(out, "% ---- FCIR RMSE: adaptive lasso | lasso | unpenalized mean | unpenalized median",
        rmse_rows(fcir,
                  list(beta_0 = "$\\bm{\\beta}_0$", B = "$\\bm{B}$", alpha_0 = "$\\bm{\\alpha}_0$", A = "$\\bm{A}$"),
                  list(list(est = "adaptive_lasso", stat = "Mean_RMSE"), list(est = "lasso", stat = "Mean_RMSE"),
                       list(est = "unpenalized", stat = "Mean_RMSE"), list(est = "unpenalized", stat = "Median_RMSE"))))

sel = read.csv(file.path(sim_dir, "FCIR", "tables", "lasso_selection_wide.csv"))
out = c(out, "", "% ---- FCIR ordinary lasso selection")
for (b in c("B", "A")) {
  d = sel[sel$Block == b, ]
  d = d[order(d$N, d$P), ]
  for (i in seq_len(nrow(d))) {
    lead = if (i == 1) sprintf("\\multirow{%d}{*}{$\\bm{%s}$}", nrow(d), b) else ""
    out = c(out, sprintf("%s & %d & %d & %d & %d & %.3f & %.3f & %.3f & %.3f \\\\", lead, d$N[i], d$P[i],
                         d$n_z[i], d$n_nz[i], d$correct_zeros[i], d$correct_nonzeros[i], d$F1[i], d$Mean_RMSE[i]))
  }
  out = c(out, "\\midrule")
}

two_col = list(list(est = "adaptive_lasso", stat = "Mean_RMSE"), list(est = "unpenalized", stat = "Mean_RMSE"),
               list(est = "unpenalized", stat = "Median_RMSE"))
fcir_i = read.csv(file.path(sim_dir, "FCIR_I", "tables", "estimator_RMSE_comparison.csv"))
out = c(out, "", "% ---- FCIR_I RMSE: adaptive lasso | unpenalized mean | unpenalized median",
        rmse_rows(fcir_i, list(Beta = "$\\mathcal{B}$", A = "$\\bm{A}$", alpha_0 = "$\\bm{\\alpha}_0$"), two_col))

fcir_m = read.csv(file.path(sim_dir, "FCIR_M", "tables", "estimator_RMSE_comparison.csv"))
out = c(out, "", "% ---- FCIR_M RMSE: adaptive lasso | unpenalized mean | unpenalized median",
        rmse_rows(fcir_m, list(beta_0 = "$\\bm{\\beta}_0$", B = "$\\bm{B}$", Theta = "$\\bm{\\Theta}$"), two_col,
                  p_first = TRUE))

# Divergence diagnostics for the text.
diag = rbind(fcir, fcir_i, fcir_m[fcir_m$Estimator == "unpenalized", ])
diag = diag[diag$Estimator == "unpenalized" & (diag$n_diverged > 0 | diag$n_nonfinite > 0), ]
out = c(out, "", "% ---- unpenalized replicates with max|est| > 10 (n_diverged) or non-finite estimates",
        sprintf("%% %s %s N=%d P=%d diverged=%d nonfinite=%d mean=%.4g median=%.4g", diag$Model, diag$Block,
                diag$N, diag$P, diag$n_diverged, diag$n_nonfinite, diag$Mean_RMSE, diag$Median_RMSE))

writeLines(out, out_file)
message("Saved: ", out_file)
