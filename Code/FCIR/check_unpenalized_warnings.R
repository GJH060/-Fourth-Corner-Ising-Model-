# Capture warnings produced by unpenalized pseudo-likelihood GLM fits.
# This script reads existing simulated data, reruns only the unpenalized
# estimation step, and saves warning logs separately. It does not overwrite
# any existing estimate Rdata files.

# ---------------- User settings ----------------
project_root <- "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
source(file.path(project_root, "Code", "FCIR", "estimate_FCIR.r"))
estimate_fun <- get("estimate_unpenalized_FCIR")

# Sparse workflow:
use_dense <- FALSE

# Start with smaller Ns if you only want a quick diagnostic, e.g. c(50, 100).
Ns <- c(50, 100, 200, 400, 800)
Ps <- c(30, 60)

# Set to a smaller number such as 100 for a quick test.
# Use NULL to check all Monte Carlo repetitions in each data file.
max_reps <- 100

output_dir <- file.path(project_root, "Simulation_Results", "FCIR", "warnings")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data_dir <- if (use_dense) {
  file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Dense")
} else {
  file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Sparse")
}

scenario_tag <- if (use_dense) "Dense" else "Sparse"

# ---------------- Warning-capture wrapper ----------------
estimate_unpenalized_warning_only <- function(Y, X, Tr) {
  warning_messages <- character(0)
  error_message <- NA_character_

  fit <- tryCatch(
    withCallingHandlers(
      estimate_fun(Y = Y, X = X, Tr = Tr),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  list(
    warnings = unique(warning_messages),
    error = error_message,
    converged = if (!is.null(fit)) fit$glm_model$converged else NA
  )
}

empty_log <- function() {
  data.frame(
    Scenario = character(),
    N = integer(),
    P = integer(),
    B = integer(),
    Type = character(),
    Message = character(),
    Converged = logical(),
    stringsAsFactors = FALSE
  )
}

# ---------------- Main loop ----------------
all_logs <- list()

for (n in Ns) {
  for (p in Ps) {
    data_file <- if (use_dense) {
      file.path(data_dir, sprintf("FCIR_data_Dense_N%d_P%d.Rdata", n, p))
    } else {
      file.path(data_dir, sprintf("FCIR_data_N%d_P%d.Rdata", n, p))
    }

    if (!file.exists(data_file)) {
      warning("Data file not found: ", data_file)
      next
    }

    message("Loading data: N=", n, ", P=", p)
    env <- new.env()
    load(data_file, envir = env)

    B_available <- dim(env$Y)[3]
    B_to_check <- if (is.null(max_reps)) B_available else min(max_reps, B_available)

    combo_logs <- vector("list", B_to_check)

    for (b in seq_len(B_to_check)) {
      if (b %% 50 == 0 || b == 1 || b == B_to_check) {
        message("Checking N=", n, ", P=", p, ", replicate ", b, "/", B_to_check)
      }

      result <- estimate_unpenalized_warning_only(
        Y = env$Y[, , b],
        X = env$X,
        Tr = env$Tr
      )

      rows <- empty_log()

      if (length(result$warnings) > 0) {
        rows <- rbind(
          rows,
          data.frame(
            Scenario = scenario_tag,
            N = n,
            P = p,
            B = b,
            Type = "warning",
            Message = result$warnings,
            Converged = result$converged,
            stringsAsFactors = FALSE
          )
        )
      }

      if (!is.na(result$error)) {
        rows <- rbind(
          rows,
          data.frame(
            Scenario = scenario_tag,
            N = n,
            P = p,
            B = b,
            Type = "error",
            Message = result$error,
            Converged = result$converged,
            stringsAsFactors = FALSE
          )
        )
      }

      combo_logs[[b]] <- rows
    }

    combo_log <- do.call(rbind, combo_logs)
    if (is.null(combo_log) || nrow(combo_log) == 0) {
      combo_log <- empty_log()
    }

    combo_base <- sprintf("unpenalized_warnings_%s_N%d_P%d", scenario_tag, n, p)
    combo_csv <- file.path(output_dir, paste0(combo_base, ".csv"))
    combo_rdata <- file.path(output_dir, paste0(combo_base, ".Rdata"))

    write.csv(combo_log, combo_csv, row.names = FALSE)
    save(combo_log, file = combo_rdata)

    all_logs[[paste0("N", n, "_P", p)]] <- combo_log
    message("Saved warning log: ", combo_csv)
  }
}

warning_log_all <- do.call(rbind, all_logs)
if (is.null(warning_log_all) || nrow(warning_log_all) == 0) {
  warning_log_all <- empty_log()
}

all_base <- paste0("unpenalized_warnings_", scenario_tag, "_ALL")
all_csv <- file.path(output_dir, paste0(all_base, ".csv"))
all_rdata <- file.path(output_dir, paste0(all_base, ".Rdata"))

write.csv(warning_log_all, all_csv, row.names = FALSE)
save(warning_log_all, file = all_rdata)

message("Finished warning capture.")
message("Combined CSV: ", all_csv)
message("Total warning/error rows: ", nrow(warning_log_all))

if (nrow(warning_log_all) > 0) {
  print(table(warning_log_all$N, warning_log_all$P))
  print(table(warning_log_all$Type, warning_log_all$Message))
}
