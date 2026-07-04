project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
rdata_dir = file.path(project_root, "Simulation_Results", "Rdata_Sparse")

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

# Helper: one-line range summary for a numeric object.
range_line <- function(name, x) {
  x <- as.numeric(x)
  n_zero <- sum(x == 0)
  sprintf("%-9s min=%8.4f  max=%8.4f  mean=%8.4f  |.|max=%8.4f  #zero=%d/%d",
          name, min(x), max(x), mean(x), max(abs(x)), n_zero, length(x))
}

# ---- The 12 fourth-corner coefficients (B_mat + A_mat) ----
# These depend only on (L, K, seed), so they are identical across all files;
# load one file and list every individual coefficient with its (l, k) index.
ref_file = file.path(rdata_dir, paste0("FCIR_data_N", Ns[1], "_P", Ps[1], ".Rdata"))
if (file.exists(ref_file)) {
  e0 = new.env(); load(ref_file, envir = e0)
  L = e0$L; K = e0$K

  coef_df = do.call(rbind, lapply(c("B", "A"), function(mat) {
    M = if (mat == "B") e0$B_mat else e0$A_mat
    data.frame(
      Matrix = mat,
      l = rep(1:L, times = K),
      k = rep(1:K, each = L),
      value = as.vector(M),
      is_zero = as.vector(M) == 0
    )
  }))

  cat("========================================================\n")
  cat("True fourth-corner coefficients (B_mat then A_mat), L=", L, " K=", K, "\n", sep = "")
  cat("--------------------------------------------------------\n")
  for (i in seq_len(nrow(coef_df))) {
    cat(sprintf("%s[%d,%d] = %8.4f%s\n",
                coef_df$Matrix[i], coef_df$l[i], coef_df$k[i],
                coef_df$value[i], ifelse(coef_df$is_zero[i], "   (zero)", "")))
  }
  cat(sprintf("\nNon-zero: %d/%d   |   Zero: %d/%d\n",
              sum(!coef_df$is_zero), nrow(coef_df),
              sum(coef_df$is_zero), nrow(coef_df)))
}

for (n in Ns) {
  for (p in Ps) {
    f = file.path(rdata_dir, paste0("FCIR_data_N", n, "_P", p, ".Rdata"))
    if (!file.exists(f)) {
      message("Missing: ", f, " - skipping.")
      next
    }
    e = new.env(); load(f, envir = e)

    cat("\n========================================================\n")
    cat(sprintf("File: FCIR_data_N%d_P%d.Rdata  (N=%d, P=%d, L=%d, K=%d, B_reps=%d)\n",
                n, p, e$N, e$P, e$L, e$K, e$B_reps))
    cat("--------------------------------------------------------\n")
    # True structural parameters (depend only on L, K, seed -> identical across files)
    cat(range_line("beta_0",  e$beta_0),  "\n")
    cat(range_line("B_mat",   e$B_mat),   "\n")
    cat(range_line("alpha_0", e$alpha_0), "\n")
    cat(range_line("A_mat",   e$A_mat),   "\n")
    # Design / trait inputs
    cat(range_line("X",  e$X),  "\n")
    cat(range_line("Tr", e$Tr), "\n")
    # Response (binary 0/1) -> report overall prevalence of 1s
    cat(sprintf("%-9s prevalence(Y==1)=%.4f  over %d cells\n",
                "Y", mean(e$Y == 1), length(e$Y)))
  }
}

cat("\nDone.\n")
