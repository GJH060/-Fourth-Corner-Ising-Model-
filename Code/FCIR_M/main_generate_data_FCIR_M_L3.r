# Regenerates the L = 3 FCIR_M simulation data used in the thesis (the
# "New_FCIR_M_..._unpenalized" tables, run before commit 8db7f4b): sparse B and
# residual edges of +/-0.25. Same generator and seed as main_generate_data_FCIR_M.r;
# the (N, P) configurations are generated in parallel because each one calls
# set.seed(seed) independently.

project_root = "D:/ising model thesis/-Fourth-Corner-Ising-Model-"
fcir_m_code_dir = file.path(project_root, "Code", "FCIR_M")
rdata_dir = file.path(project_root, "Simulation_Results", "FCIR_M", "Rdata")

Ns = c(50, 100, 200, 400, 800)
Ps = c(10, 20)

L = 3
K = 2
B_reps = 1000
seed = 2026
sparse_B = TRUE
edge_magnitude = 0.25

if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

grid = expand.grid(N = Ns, P = Ps)
grid$file = file.path(rdata_dir, paste0("FCIR_M_data_L", L, "_N", grid$N, "_P", grid$P, ".Rdata"))
grid = grid[!file.exists(grid$file), ]
# Largest configurations first so they do not finish last on their own.
grid = grid[order(-grid$N * grid$P), ]

total_start = Sys.time()

if (nrow(grid) > 0) {
  cl = parallel::makeCluster(nrow(grid))
  parallel::clusterCall(cl, function(dir) {
    source(file.path(dir, "generate_FCIR_M.r"))
  }, dir = fcir_m_code_dir)

  parallel::clusterApplyLB(cl, split(grid, seq_len(nrow(grid))),
                           function(g, L, K, B_reps, seed, sparse_B, edge_magnitude) {
    generate_fcir_M_data(N = g$N, P = g$P, L = L, K = K, B_reps = B_reps,
                         seed = seed, filename = g$file,
                         sparse_B = sparse_B, edge_magnitude = edge_magnitude)
    g$file
  }, L = L, K = K, B_reps = B_reps, seed = seed,
     sparse_B = sparse_B, edge_magnitude = edge_magnitude)

  parallel::stopCluster(cl)
}

print(paste("Total FCIR_M (L = 3) data generation time:", format(Sys.time() - total_start)))
