project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_i_code_dir = file.path(project_root, "Code", "FCIR_I")

# Set to TRUE to fit the sparse-Beta data, or FALSE for the original dense
# Beta_mat data. Must match the setting used in main_generate_data_FCIR_I.r.
use_sparse_beta = TRUE

if (use_sparse_beta) {
  rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata_Sparse_Beta")
  est_tag = "unpenalized_sparse_Beta"
  data_label = "FCIR_I sparse-Beta"
  beta_structure = "sparse"
} else {
  rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata")
  est_tag = "unpenalized"
  data_label = "FCIR_I"
  beta_structure = "dense"
}

source(file.path(fcir_i_code_dir, "estimate_FCIR_I.r"))

if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Package 'callr' is required. Install it with install.packages('callr').")
}

# This machine suffers sporadic native crashes (access violations inside R.dll)
# that kill the whole R process, so each batch runs in a throwaway child
# process: a crash takes down only that child and the batch is retried.
# Finished batches are checkpointed, so an interrupted run resumes where it
# stopped instead of restarting from the first replicate.
n_cores = min(6, max(1, parallel::detectCores() - 1))
batch_size = 50
max_attempts = 5

checkpoint_root = file.path(rdata_dir, "checkpoints")

fit_batch = function(data_file, batch_idx, code_dir, n_cores) {
  source(file.path(code_dir, "estimate_FCIR_I.r"))

  # Keep only this batch's slices plus X/Tr, so the full Y array is never
  # exported to the workers.
  data_env = new.env()
  load(data_file, envir = data_env)
  X = data_env$X
  Tr = data_env$Tr
  Y_batch = asplit(data_env$Y[, , batch_idx, drop = FALSE], 3)
  rm(data_env)

  fit_one = function(Y_b) {
    res = suppressWarnings(
      get("estimate_unpenalized_FCIR_I")(Y = Y_b, X = X, Tr = Tr)
    )
    list(Beta_mat = res$Beta_mat,
         alpha_0 = res$alpha_0,
         A_mat = res$A_mat,
         converged = res$glm_model$converged)
  }

  if (n_cores <= 1) {
    return(lapply(Y_batch, fit_one))
  }

  library(foreach)
  library(doParallel)
  cl = parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  registerDoParallel(cl)
  parallel::clusterCall(cl, function(dir) {
    source(file.path(dir, "estimate_FCIR_I.r"))
  }, dir = code_dir)

  foreach(Y_b = Y_batch, .export = c("X", "Tr")) %dopar% fit_one(Y_b)
}

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(
      rdata_dir,
      paste0("Rescale_FCIR_I_data_N", n, "_P", p, ".Rdata")
    )
    legacy_data_filename = file.path(
      rdata_dir,
      paste0("FCIR_I_data_N", n, "_P", p, ".Rdata")
    )
    if (!file.exists(data_filename) && file.exists(legacy_data_filename)) {
      data_filename = legacy_data_filename
    }
    est_prefix = if (identical(data_filename, legacy_data_filename)) "" else "Rescale_"
    est_filename = file.path(
      rdata_dir,
      paste0(est_prefix, "FCIR_I_estimates_", est_tag,
             "_N", n, "_P", p, ".Rdata")
    )

    if (!file.exists(data_filename)) {
      print(paste(data_label, "data not found, skipping:", data_filename))
      next
    }

    if (file.exists(est_filename)) {
      print(paste(data_label, "unpenalized estimates already exist for N =", n,
                  ", P =", p, "- Skipping."))
      next
    }

    print(paste("Loading", data_label, "data ( N =", n, ", P =", p, ")..."))
    load(data_filename)   # Y, X, Tr, Beta_mat, alpha_0, A_mat, ...

    class_counts = t(vapply(seq_len(dim(Y)[3]), function(b) {
      c(zeros = sum(Y[, , b] == 0), ones = sum(Y[, , b] == 1))
    }, numeric(2)))
    bad_reps = which(class_counts[, "zeros"] == 0 | class_counts[, "ones"] == 0)
    if (length(bad_reps) > 0) {
      stop(paste0(
        data_label, " data contain ", length(bad_reps),
        " replicate(s) with only one response class at N=", n, ", P=", p,
        ". First bad reps: ", paste(head(bad_reps, 20), collapse = ", "),
        ". Regenerate the data before fitting."
      ))
    }

    print(paste("Starting", B_reps, data_label, "unpenalized fits for N =",
                n, "P =", p, "..."))

    est_Beta_mat = array(NA, dim = c(P, L, B_reps))
    est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
    est_A_mat = array(NA, dim = c(L, K, B_reps))
    est_converged = logical(B_reps)

    rm(Y)

    checkpoint_dir = file.path(checkpoint_root,
                               sub("\\.Rdata$", "", basename(est_filename)))
    if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

    for (batch_start in seq(1, B_reps, by = batch_size)) {
      batch_idx = batch_start:min(batch_start + batch_size - 1, B_reps)
      checkpoint_file = file.path(checkpoint_dir,
                                  paste0("batch_", batch_start, ".rds"))

      if (file.exists(checkpoint_file)) {
        batch_results = readRDS(checkpoint_file)
      } else {
        batch_results = NULL
        for (attempt in seq_len(max_attempts)) {
          batch_results = tryCatch(
            callr::r(fit_batch,
                     args = list(data_file = data_filename,
                                 batch_idx = batch_idx,
                                 code_dir = fcir_i_code_dir,
                                 n_cores = n_cores)),
            error = function(e) {
              print(paste0("Batch ", batch_start, " attempt ", attempt, "/",
                           max_attempts, " failed: ", conditionMessage(e)))
              NULL
            }
          )
          if (!is.null(batch_results)) break
        }

        if (is.null(batch_results)) {
          stop(paste0("Batch starting at replicate ", batch_start,
                      " failed ", max_attempts, " times at N=", n, ", P=", p,
                      ". Completed batches are checkpointed in ", checkpoint_dir,
                      "; re-run this script to resume."))
        }
        saveRDS(batch_results, checkpoint_file)
      }

      for (i in seq_along(batch_idx)) {
        b = batch_idx[i]
        est_Beta_mat[, , b] = batch_results[[i]]$Beta_mat
        est_alpha_0[b, ] = batch_results[[i]]$alpha_0
        est_A_mat[, , b] = batch_results[[i]]$A_mat
        est_converged[b] = batch_results[[i]]$converged
      }
      rm(batch_results)
      print(paste("Completed replicates", max(batch_idx), "of", B_reps))
    }

    method = "unpenalized"
    save(est_Beta_mat, est_alpha_0, est_A_mat, est_converged,
         Beta_mat, alpha_0, A_mat,
         N = n, P = p, L, K, B_reps, seed, method, beta_structure,
         file = est_filename)
    print(paste("Saved:", est_filename,
                "- converged:", sum(est_converged), "/", B_reps))
    unlink(checkpoint_dir, recursive = TRUE)
  }
}

print(paste("Total", data_label, "unpenalized execution time:",
            Sys.time() - total_start))
