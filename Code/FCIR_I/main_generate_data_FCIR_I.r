project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_i_code_dir = file.path(project_root, "Code", "FCIR_I")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_I", "Rdata")
source(file.path(fcir_i_code_dir, "generate_FCIR_I.r"))

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
      print(paste("Generating FCIR_I data ( N =", n, ", P =", p, ")..."))
      generate_fcir_I_data(N = n, P = p, L = L, K = K, B_reps = B_reps,
                           seed = seed, filename = data_filename)
      print(paste("Saved:", data_filename))
    } else {
      print(paste("FCIR_I data already exists for N =", n, ", P =", p, "- Skipping generation."))
    }
  }
}

print(paste("Total FCIR_I data generation time:", Sys.time() - total_start))
