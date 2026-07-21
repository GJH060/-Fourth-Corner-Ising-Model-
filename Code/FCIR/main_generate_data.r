project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_code_dir = file.path(project_root, "Code", "FCIR")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Sparse")
source(file.path(fcir_code_dir, "generate_FCIR.r"))

# 1. Define the parameter grid you want to sweep across
Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

# Global fixed parameters
L = 3
K = 2
B_reps = 1000
seed = 2026

if (!dir.exists(rdata_dir)) {
  dir.create(rdata_dir, recursive = TRUE)
}

total_start = Sys.time()

# 2. Generate data only. Estimation is handled by main.r.
for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(rdata_dir, paste0("FCIR_data_N", n, "_P", p, ".Rdata"))

    if (!file.exists(data_filename)) {
      print(paste("Generating data ( N =", n, ", P =", p, ")..."))
      generate_dense_fcir_data(N = n, P = p, L = L, K = K, B_reps = B_reps,
                               seed = seed, filename = data_filename)
      print(paste("Saved:", data_filename))
    } else {
      print(paste("Data already exists for N =", n, ", P =", p, "- Skipping generation."))
    }
  }
}

print(paste("Total data generation time:", Sys.time() - total_start))
