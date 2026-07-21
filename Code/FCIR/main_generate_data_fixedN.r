project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_code_dir = file.path(project_root, "Code", "FCIR")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Sparse")
source(file.path(fcir_code_dir, "generate_FCIR.r"))

# 1. Fixed sample size, sweep the number of species P.
N_fixed = 200
Ps = c(30, 60, 120, 240, 480)

# Global fixed parameters
L = 3
K = 2
B_reps = 1000
seed = 42

if (!dir.exists(rdata_dir)) {
  dir.create(rdata_dir, recursive = TRUE)
}

total_start = Sys.time()

# 2. Generate data only (fixed N = 200). Estimation is handled by main.r.
for (p in Ps) {
  data_filename = file.path(rdata_dir, paste0("FCIR_data_N", N_fixed, "_P", p, ".Rdata"))

  if (!file.exists(data_filename)) {
    print(paste("Generating data ( N =", N_fixed, ", P =", p, ")..."))
    generate_dense_fcir_data(N = N_fixed, P = p, L = L, K = K, B_reps = B_reps,
                             seed = seed, filename = data_filename)
    print(paste("Saved:", data_filename))
  } else {
    print(paste("Data already exists for N =", N_fixed, ", P =", p, "- Skipping generation."))
  }
}

print(paste("Total data generation time:", Sys.time() - total_start))
