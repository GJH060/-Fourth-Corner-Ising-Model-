# Render every Application result CSV as a bordered ASCII table, in the same
# style as Simulation_Results/FCIR/tables/*.txt.
#
#   Rscript Code/Application/format_app_tables.r
#   source("Code/Application/format_app_tables.r")
#
# Writes <name>.txt beside each <name>.csv. The CSVs are left untouched.

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
table_dir = file.path(project_root, "Simulation_Results", "Application", "tables")

# Verbatim from Code/FCIR/evaluate_selection_adaptive_lasso.r so the borders,
# padding and alignment match the simulation tables exactly.
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

# The simulation tables all round to 3 decimals, which is fine when every number
# sits near 1. Application tables mix scales -- lambda_min reaches 1.2e-4 and
# would round to a bare "0" at 3 decimals, while pct_selected is 74.7. Pick the
# decimals per column so the smallest non-zero value in it still shows `sig`
# significant digits.
choose_digits <- function(col, sig = 3, max_digits = 6) {
  v <- abs(col[is.finite(col) & col != 0])
  if (length(v) == 0) return(0)
  if (all(v == round(v))) return(0)                 # integer-valued column
  d <- ceiling(-log10(min(v))) + sig - 1
  min(max(d, 1), max_digits)
}

# format_ascii_table() takes one `digits` for the whole frame, so round each
# numeric column to its own precision first and hand over pre-rounded values.
round_per_column <- function(df, sig = 3, max_digits = 6) {
  for (nm in names(df)) {
    if (is.numeric(df[[nm]])) {
      df[[nm]] <- round(df[[nm]], choose_digits(df[[nm]], sig, max_digits))
    }
  }
  df
}

csv_files = sort(list.files(table_dir, pattern = "[.]csv$", full.names = TRUE))
if (length(csv_files) == 0) stop("no CSVs found in ", table_dir)

for (f in csv_files) {
  df = read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                fileEncoding = "UTF-8-BOM")
  lines = format_ascii_table(round_per_column(df), digits = 6)
  out = sub("[.]csv$", ".txt", f)
  writeLines(lines, out)
  cat(sprintf("  %-38s %4d rows x %2d cols  ->  %s (%d chars wide)\n",
              basename(f), nrow(df), ncol(df), basename(out), nchar(lines[1])))
}

cat("\nWrote ", length(csv_files), " tables to ", table_dir, "\n", sep = "")
