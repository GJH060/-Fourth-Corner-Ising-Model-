# Scans every saved .Rdata under Simulation_Results for blatant signs of the
# memory corruption traced to this machine (NaN/Inf entries, absurd magnitudes,
# arrays whose dimensions disagree with the N/P/L/K/B_reps stored beside them).
#
# A flipped mantissa bit produces a perfectly plausible number, so a clean
# report here does NOT prove the results are sound -- only a re-run on healthy
# hardware can do that. This catches the corruption severe enough to show.

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
results_root = file.path(project_root, "Simulation_Results")
report_file = file.path(results_root, "result_integrity_scan.csv")

# Logistic fits under separation can legitimately produce large coefficients,
# but nothing in this model family approaches 1e100.
magnitude_limit = 1e100

summarise_numeric <- function(v) {
  v = as.numeric(v)
  finite = is.finite(v)
  list(
    n = length(v),
    n_na = sum(is.na(v)) - sum(is.nan(v)),
    n_nan = sum(is.nan(v)),
    n_inf = sum(is.infinite(v)),
    max_abs = if (any(finite)) max(abs(v[finite])) else NA_real_
  )
}

check_dims <- function(e) {
  problems = character(0)
  needed = c("P", "L", "K", "B_reps")
  has_meta = vapply(needed, exists, logical(1), envir = e, inherits = FALSE)
  if (!all(has_meta)) return(problems)

  expect <- function(nm, want) {
    if (!exists(nm, envir = e, inherits = FALSE)) return(invisible(NULL))
    obj = get(nm, envir = e)
    got = if (is.null(dim(obj))) length(obj) else dim(obj)
    if (!identical(as.integer(got), as.integer(want))) {
      problems <<- c(problems, paste0(nm, " is ", paste(got, collapse = "x"),
                                      " but should be ", paste(want, collapse = "x")))
    }
  }

  P = e$P; L = e$L; K = e$K; B = e$B_reps
  expect("est_Beta_mat", c(P, L, B))
  expect("est_alpha_0", c(B, L))
  expect("est_A_mat", c(L, K, B))
  expect("est_converged", B)
  expect("Beta_mat", c(P, L))
  expect("A_mat", c(L, K))
  problems
}

files = list.files(results_root, pattern = "\\.Rdata$",
                   recursive = TRUE, full.names = TRUE)
if (length(files) == 0) {
  stop("No .Rdata files found under ", results_root)
}

message("Scanning ", length(files), " .Rdata files...")

rows = vector("list", length(files))

for (i in seq_along(files)) {
  f = files[i]
  rel = substring(f, nchar(results_root) + 2)

  e = new.env()
  load_error = tryCatch({
    load(f, envir = e)
    NA_character_
  }, error = function(err) conditionMessage(err))

  if (!is.na(load_error)) {
    rows[[i]] = data.frame(file = rel, n_values = NA_integer_,
                           n_na = NA_integer_, n_nan = NA_integer_,
                           n_inf = NA_integer_, max_abs = NA_real_,
                           suspect_objects = "", dim_problems = "",
                           load_error = load_error,
                           stringsAsFactors = FALSE)
    message(sprintf("[%d/%d] %s -- FAILED TO LOAD", i, length(files), rel))
    next
  }

  n_values = 0; n_na = 0; n_nan = 0; n_inf = 0
  max_abs = NA_real_
  suspect = character(0)

  for (nm in ls(e)) {
    obj = get(nm, envir = e)
    if (!is.numeric(obj)) next
    s = summarise_numeric(obj)
    n_values = n_values + s$n
    n_na = n_na + s$n_na
    n_nan = n_nan + s$n_nan
    n_inf = n_inf + s$n_inf
    max_abs = max(max_abs, s$max_abs, na.rm = TRUE)
    if (s$n_nan > 0 || s$n_inf > 0 ||
        (!is.na(s$max_abs) && s$max_abs > magnitude_limit)) {
      suspect = c(suspect, nm)
    }
  }

  dim_problems = check_dims(e)

  rows[[i]] = data.frame(
    file = rel,
    n_values = n_values,
    n_na = n_na,
    n_nan = n_nan,
    n_inf = n_inf,
    max_abs = max_abs,
    suspect_objects = paste(suspect, collapse = "; "),
    dim_problems = paste(dim_problems, collapse = "; "),
    load_error = NA_character_,
    stringsAsFactors = FALSE
  )

  message(sprintf("[%d/%d] %s", i, length(files), rel))
  rm(e)
  gc(verbose = FALSE)
}

report = do.call(rbind, rows)
write.csv(report, report_file, row.names = FALSE)

flagged = report[!is.na(report$load_error) |
                   report$suspect_objects != "" |
                   report$dim_problems != "", ]

message("")
message("Scanned ", nrow(report), " files. Full report: ", report_file)

if (nrow(flagged) == 0) {
  message("No NaN/Inf entries, no absurd magnitudes, no dimension mismatches.")
  message("Note: this cannot rule out silently altered values -- re-run on ",
          "healthy hardware and compare to confirm.")
} else {
  message("FLAGGED ", nrow(flagged), " file(s):")
  print(flagged[, c("file", "n_nan", "n_inf", "max_abs",
                    "suspect_objects", "dim_problems", "load_error")],
        row.names = FALSE)
}

# Plain NA is expected: glm returns NA coefficients for aliased columns.
message("")
message("Files containing NA coefficients (usually legitimate aliasing): ",
        sum(report$n_na > 0, na.rm = TRUE))
