project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_i_code_dir = file.path(project_root, "Code", "FCIR_I")

# Set to TRUE to generate data with a sparse Beta_mat, or FALSE for the
# original dense Beta_mat.
use_sparse_beta = TRUE

if (use_sparse_beta) {
  rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata_Sparse_Beta")
  source(file.path(fcir_i_code_dir, "generate_sparse_Beta_FCIR_I.r"))
  generate_data = generate_sparse_beta_fcir_I_data
  data_label = "FCIR_I sparse-Beta"
} else {
  rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata")
  source(file.path(fcir_i_code_dir, "generate_FCIR_I.r"))
  generate_data = generate_fcir_I_data
  data_label = "FCIR_I"
}

Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

L = 3
K = 2
B_reps = 1000
seed = 2026

if (!dir.exists(rdata_dir)) {
  dir.create(rdata_dir, recursive = TRUE)
}

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("FCIR_I_data_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("Generating", data_label, "data ( N =", n, ", P =", p, ")..."))
      generate_data(N = n, P = p, L = L, K = K, B_reps = B_reps,
                    seed = seed, filename = data_filename)
      print(paste("Saved:", data_filename))
    } else {
      print(paste(data_label, "data already exists for N =", n, ", P =", p, "- Skipping generation."))
    }
  }
}

print(paste("Total", data_label, "data generation time:", Sys.time() - total_start))
