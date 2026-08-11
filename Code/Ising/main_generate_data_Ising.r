source("F:/ising model thesis/-Fourth-Corner-Ising-Model-/Code/Ising/ising_config.r")
source(file.path(ising_code_dir, "generate_Ising.r"))

if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

total_start = Sys.time()

for (n in Ns) {
  for (p in Ps) {
    data_filename = file.path(
      rdata_dir, paste0("Ising_data", setting_tag, "_N", n, "_P", p, ".Rdata")
    )

    if (file.exists(data_filename)) {
      print(paste("Ising data already exists for N =", n, ", P =", p,
                  "- Skipping generation."))
      next
    }

    print(paste("Generating Ising data ( N =", n, ", P =", p, ")..."))
    generate_ising_data(N = n, P = p, B_reps = B_reps, seed = seed,
                        filename = data_filename,
                        theta_jj_range = theta_jj_range,
                        theta_jj_density = theta_jj_density,
                        edge_density = edge_density,
                        edge_values = edge_values,
                        edge_magnitude = edge_magnitude,
                        theta_jjp_range = theta_jjp_range,
                        nIter = nIter)
    print(paste("Saved:", data_filename))
  }
}

print(paste("Total Ising data generation time:", Sys.time() - total_start))
