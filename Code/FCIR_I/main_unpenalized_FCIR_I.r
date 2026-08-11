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

library(foreach)
library(doParallel)

n_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(n_cores)
registerDoParallel(cl)
clusterCall(cl, function(dir) {
  source(file.path(dir, "estimate_FCIR_I.r"))
}, dir = fcir_i_code_dir)

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("FCIR_I_data_N", n, "_P", p, ".Rdata"))
    est_filename = file.path(rdata_dir, paste0("FCIR_I_estimates_", est_tag, "_N", n, "_P", p, ".Rdata"))

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

    Y_slices = asplit(Y, 3)
    rm(Y)

    unpenalized_results = foreach(Y_b = Y_slices) %dopar% {
      res = estimate_unpenalized_FCIR_I(Y = Y_b, X = X, Tr = Tr)
      list(Beta_mat = res$Beta_mat,
           alpha_0 = res$alpha_0,
           A_mat = res$A_mat,
           converged = res$glm_model$converged)
    }
    rm(Y_slices)

    for (b in 1:B_reps) {
      est_Beta_mat[, , b] = unpenalized_results[[b]]$Beta_mat
      est_alpha_0[b, ] = unpenalized_results[[b]]$alpha_0
      est_A_mat[, , b] = unpenalized_results[[b]]$A_mat
      est_converged[b] = unpenalized_results[[b]]$converged
    }
    rm(unpenalized_results)

    method = "unpenalized"
    save(est_Beta_mat, est_alpha_0, est_A_mat, est_converged,
         Beta_mat, alpha_0, A_mat,
         N = n, P = p, L, K, B_reps, seed, method, beta_structure,
         file = est_filename)
    print(paste("Saved:", est_filename,
                "- converged:", sum(est_converged), "/", B_reps))
  }
}

stopCluster(cl)
print(paste("Total", data_label, "unpenalized execution time:",
            Sys.time() - total_start))
