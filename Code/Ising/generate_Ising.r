library(IsingSampler)

generate_ising_data <- function(N, P, B_reps, seed, filename,
                                theta_jj_range = c(-0.5, 0.5),
                                theta_jj_density = 1 / 3,
                                edge_density = 1 / 3,
                                edge_values = "two_point",
                                edge_magnitude = 0.25,
                                theta_jjp_range = c(-0.5, 0.5),
                                nIter = 1000) {
  # Plain Ising model: no environment matrix, no trait matrix. Parameters are
  # the node thresholds theta_jj and the pairwise interactions theta_jj', both
  # drawn once and held fixed across the B_reps datasets.
  #
  # Thresholds: a fraction theta_jj_density of the P nodes are non-zero
  # (Unif(theta_jj_range)); the rest are exactly zero. Density 1 recovers the
  # earlier dense-threshold design.
  #
  # The interaction network follows the FCIR_M convention (generate_FCIR_M.r):
  # a fraction edge_density of the P*(P-1)/2 edges is non-zero, each equal to
  # +/- edge_magnitude. Set edge_values = "uniform" to draw the non-zero edges
  # from Unif(theta_jjp_range) instead, which gives a continuous spread of edge
  # strengths including weak edges that are much harder to select.

  if (!edge_values %in% c("two_point", "uniform")) {
    stop("edge_values must be 'two_point' or 'uniform'.")
  }
  if (theta_jj_density < 0 || theta_jj_density > 1) {
    stop("theta_jj_density must be in [0, 1].")
  }

  set.seed(seed)

  # Thresholds: sparse or dense according to theta_jj_density.
  theta_jj = numeric(P)
  n_nonzero_jj = round(P * theta_jj_density)
  if (n_nonzero_jj > 0) {
    sel_jj = sample.int(P, size = n_nonzero_jj, replace = FALSE)
    theta_jj[sel_jj] = runif(n_nonzero_jj, theta_jj_range[1], theta_jj_range[2])
  }

  # Pairwise interactions: symmetric with a zero diagonal.
  Theta = matrix(0, nrow = P, ncol = P)
  n_edges = P * (P - 1) / 2
  # round(), not floor(): n_edges * (1/3) can land just below the integer in
  # floating point (e.g. 45 * (1/3) = 14.999...), which would silently drop an
  # edge relative to the floor(n_edges / 3) form used in generate_FCIR_M.r.
  n_nonzero = round(n_edges * edge_density)
  sel_nonzeros = sample.int(n_edges, size = n_nonzero, replace = FALSE)

  edge_vec = numeric(n_edges)
  edge_vec[sel_nonzeros] = if (edge_values == "two_point") {
    sample(c(-edge_magnitude, edge_magnitude), size = n_nonzero, replace = TRUE)
  } else {
    runif(n_nonzero, theta_jjp_range[1], theta_jjp_range[2])
  }

  Theta[upper.tri(Theta)] = edge_vec
  Theta[lower.tri(Theta)] = t(Theta)[lower.tri(Theta)]

  # Without site covariates every site is an i.i.d. draw from the same Ising
  # distribution, so one IsingSampler call yields the whole N x P matrix.
  Y = array(NA, dim = c(N, P, B_reps))
  for (b in 1:B_reps) {
    Y[, , b] = IsingSampler(N, Theta, theta_jj, beta = 1, nIter = nIter,
                            responses = c(0L, 1L), method = "MH")
  }

  save(Y, theta_jj, Theta, N, P, B_reps, seed,
       theta_jj_range, theta_jj_density, edge_density, edge_values,
       edge_magnitude, theta_jjp_range, nIter,
       file = filename)
}
