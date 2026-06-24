library(glmnet)
library(dplyr)

project_root  = "E:/RStudio/thesis"
fcir_code_dir = file.path(project_root, "Code", "FCIR")
source(file.path(fcir_code_dir, "estimate_adaptive_lasso_FCIR.r"))

# ---- Settings (edit as needed) --------------------------------------------
use_dense   = FALSE                 # FALSE -> Rdata_Sparse, TRUE -> Rdata_Dense
init_method = "unpenalized"         # must be a valid init for estimate_adaptive_lasso_FCIR
gamma_value = 1
trace_cases = list(                 # which (N, P) datasets to inspect
  c(N = 50,  P = 30),
  c(N = 200, P = 30),
  c(N = 800, P = 30)
)
trace_reps  = 1:5                    # which replicate indices b to trace
tol         = 1e-8                   # |true coef| <= tol counts as a true zero

rdata_dir = file.path(project_root, "Simulation_Results",
                      if (use_dense) "Rdata_Dense" else "Rdata_Sparse")
infix     = if (use_dense) "_Dense" else ""
out_dir   = file.path(project_root, "Simulation_Results", "tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Helper: render a data.frame as a bordered ASCII table -----------------
format_ascii_table <- function(x, digits = 3) {
  x <- as.data.frame(x)
  formatted <- lapply(x, function(col) {
    if (is.numeric(col)) {
      out <- format(round(col, digits), trim = TRUE, scientific = FALSE)
    } else {
      out <- as.character(col)
    }
    out[is.na(out)] <- ""
    out
  })
  formatted <- as.data.frame(formatted, stringsAsFactors = FALSE, check.names = FALSE)
  names(formatted) <- names(x)

  widths <- pmax(
    nchar(names(formatted)),
    vapply(formatted, function(col) max(nchar(col)), integer(1))
  )
  border <- paste0("+", paste(vapply(widths, function(w) paste(rep("-", w + 2), collapse = ""), character(1)), collapse = "+"), "+")
  row_line <- function(vals) {
    paste0("| ", paste(mapply(function(val, w) sprintf(paste0("%-", w, "s"), val), vals, widths), collapse = " | "), " |")
  }

  c(
    border,
    row_line(names(formatted)),
    border,
    apply(formatted, 1, row_line),
    border
  )
}

build_block_table <- function(beta_init, pen_factor, beta_0, B_mat, alpha_0, A_mat) {
  L <- length(beta_0); K <- ncol(B_mat)

  block  <- c(rep("beta_0",  L - 1),
              rep("B",       L * K),
              rep("alpha_0", L),
              rep("A",       L * K))
  # element label within each block
  elem   <- c(paste0("beta_0[", 2:L, "]"),
              paste0("B[",  rep(1:L, times = K), ",", rep(1:K, each = L), "]"),
              paste0("alpha_0[", 1:L, "]"),
              paste0("A[",  rep(1:L, times = K), ",", rep(1:K, each = L), "]"))
  true_v <- c(beta_0[2:L],
              as.vector(B_mat),
              alpha_0,
              as.vector(A_mat))

  data.frame(
    block       = block,
    element     = elem,
    true_value  = true_v,
    true_zero   = abs(true_v) <= tol,
    beta_init   = beta_init,
    pen_factor  = pen_factor,
    stringsAsFactors = FALSE
  )
}

# ---- Main tracing loop -----------------------------------------------------
all_rows = list()

for (case in trace_cases) {
  n = as.integer(case[["N"]]); p = as.integer(case[["P"]])
  data_filename = file.path(rdata_dir, paste0("FCIR_data", infix, "_N", n, "_P", p, ".Rdata"))

  if (!file.exists(data_filename)) {
    warning("Data not found, skipping: ", data_filename)
    next
  }

  e = new.env(); load(data_filename, envir = e)

  for (b in trace_reps) {
    if (b > dim(e$Y)[3]) next

    res = estimate_adaptive_lasso_FCIR(
      Y = e$Y[, , b], X = e$X, Tr = e$Tr,
      gamma = gamma_value, init = init_method,
      lambda = "lambda.min", use_cv = TRUE, cv_group_by_site = TRUE
    )

    tab = build_block_table(res$beta_init, res$penalty_factor,
                            e$beta_0, e$B_mat, e$alpha_0, e$A_mat)
    tab$N = n; tab$P = p; tab$rep = b
    all_rows[[paste(n, p, b, sep = "_")]] = tab
  }
  message("Traced N=", n, ", P=", p)
}

if (length(all_rows) == 0) stop("No datasets traced. Check rdata_dir / trace_cases.")

trace_df = bind_rows(all_rows)

# ---- Per-element raw trace (saved) ----------------------------------------
raw_csv = file.path(out_dir, paste0("adaptive_lasso_weight_trace_", init_method,
                                    if (use_dense) "_Dense" else "_Sparse", ".csv"))
write.csv(trace_df, raw_csv, row.names = FALSE)
message("Saved per-element trace: ", raw_csv)

# ---- Block-level summary, split by true zero vs non-zero -------------------
# This is the key diagnostic view: for each (N, P, block, true_zero) bucket,
# average beta_init magnitude and pen_factor across the traced replicates.
block_summary = trace_df %>%
  group_by(N, P, block, true_zero) %>%
  summarise(
    n_elem            = n(),
    mean_abs_beta_init = mean(abs(beta_init)),
    mean_pen_factor    = mean(pen_factor),
    median_pen_factor  = median(pen_factor),
    .groups = "drop"
  ) %>%
  arrange(N, P, block, true_zero)

scenario_tag = if (use_dense) "_Dense" else "_Sparse"
summary_csv = file.path(out_dir, paste0("adaptive_lasso_weight_trace_summary_", init_method,
                                        scenario_tag, ".csv"))
write.csv(block_summary, summary_csv, row.names = FALSE)
message("Saved block-level summary: ", summary_csv)

summary_txt = file.path(out_dir, paste0("adaptive_lasso_weight_trace_summary_", init_method,
                                        scenario_tag, ".txt"))
writeLines(format_ascii_table(block_summary), summary_txt)
message("Saved bordered block-level summary: ", summary_txt)

writeLines(format_ascii_table(block_summary))
message("Finished adaptive lasso weight tracing.")
