# Shared configuration for the plain-Ising pipeline. All three main_* scripts
# source this file so the grid and the parameter settings cannot drift apart
# between generation, estimation and evaluation.

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
ising_code_dir = file.path(project_root, "Code", "Ising_symmetrizingapproach_NOTUSED")
rdata_dir = file.path(project_root, "Simulation_Results", "Ising", "Rdata")
table_dir = file.path(project_root, "Simulation_Results", "Ising", "tables")

Ns = c(50, 100, 200, 400, 800)
Ps = c(10, 20)

B_reps = 200
seed = 2026

# --- Ising parameters -------------------------------------------------------
# Thresholds: a fraction theta_jj_density of nodes are non-zero, each drawn
# from Unif(theta_jj_range); the rest are exactly zero. Set density = 1 for
# the earlier dense-threshold design.
theta_jj_range = c(-0.5, 0.5)
theta_jj_density = 1 / 3

# Pairwise interactions follow the FCIR_M convention (generate_FCIR_M.r): one
# third of the P*(P-1)/2 edges are non-zero, each drawn from {-0.25, +0.25}.
edge_density = 1/3
edge_values = "two_point"        # "two_point" or "uniform"
edge_magnitude = 0.5          # used when edge_values = "two_point"
theta_jjp_range = c(-0.5, 0.5)   # used when edge_values = "uniform"

nIter = 1000

# --- Adaptive lasso settings ------------------------------------------------
gamma_value = 1
# "ridge" rather than "unpenalized": the node-wise unpenalized MLE is separated
# whenever N is small relative to P, which drives the adaptive weights to zero.
init_method = "unpenalized"
lambda_rule = "lambda.min"
symmetrize_rule = "and"

# Filename tag: empty only for the original default (dense thresholds +
# edge_density = 1/3, two_point). Sparse-threshold runs get a jjdens* suffix so
# they sit alongside the old Rdata instead of overwriting it.
ising_setting_tag <- function(edge_density, edge_values, theta_jj_density = 1) {
  dens_str <- function(d) sub("\\.", "", format(round(d, 3), nsmall = 3))
  parts = character(0)
  if (!(isTRUE(all.equal(edge_density, 1 / 3)) && edge_values == "two_point")) {
    parts = c(parts, paste0(edge_values, "_dens", dens_str(edge_density)))
  }
  if (!isTRUE(all.equal(theta_jj_density, 1))) {
    parts = c(parts, paste0("jjdens", dens_str(theta_jj_density)))
  }
  if (length(parts) == 0) return("")
  paste0("_", paste(parts, collapse = "_"))
}

setting_tag = ising_setting_tag(edge_density, edge_values, theta_jj_density)
