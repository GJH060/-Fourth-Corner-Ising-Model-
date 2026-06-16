# Fourth Corner Ising Regression Model

Simulation code for the **Fourth-Corner Ising Regression (FCIR)** model — a
site-varying Ising model whose main effects `theta_{s,jj}` and pairwise
interactions `theta_{s,jj'}` are driven by environmental covariates `x_s` and
species traits `t_j` (see `report.tex` for the full methodology). Models are
fit by **maximum pseudo-likelihood** (a large logistic regression via
`glmnet`).

## Repository layout

```
Code/FCIR/              R source (data generation, estimation, analysis, plots)
Simulation_Results/
  Rdata_Sparse/         simulated data + estimates, sparse-truth scenario   (git-ignored)
  Rdata_Dense/          simulated data + estimates, fully-dense scenario     (git-ignored)
  plots/                all generated figures
  warnings/             GLM fit warning logs                                 (git-ignored)
report.tex              manuscript (model definition and notation)
```

`.Rdata` outputs are large binaries and are **git-ignored**; regenerate them by
running the `main` scripts below.

## Code/FCIR scripts

### Data generation
| File | Purpose |
|------|---------|
| `generate_FCIR.r` | Generate **sparse-truth** FCIR data via Ising MCMC (`IsingSampler`). `B`/`A` have a fraction of true zeros (`prob_zero`). |
| `generate_without_sparsity.r` | Generate **fully-dense** FCIR data — all of `beta_0, B, alpha_0, A` are non-zero. |

### Estimation
| File | Purpose |
|------|---------|
| `estimate_FCIR.r` | `estimate_unpenalized_FCIR()` — unpenalized MPLE (plain logistic regression on the pseudo-likelihood design). |
| `estimate_penalized_FCIR.r` | `estimate_penalized_FCIR()` — penalized MPLE via `glmnet` (ridge/lasso/elastic-net; fixed `lambda` or `cv.glmnet`; supports **site-grouped CV folds** so all `P` rows of a site stay in one fold). |

### Simulation drivers
| File | Purpose |
|------|---------|
| `main.r` | **Sparse** scenario sweep over `N x P`. For each cell: generate data, run unpenalized MPLE, fixed-`lambda` ridge (with `lambda` scaled by `1/sqrt(N*P)`), and site-grouped CV ridge. Writes to `Rdata_Sparse/`. |
| `main_without_sparsity.r` | Same pipeline for the **dense** scenario. Writes to `Rdata_Dense/`. |

### Analysis & convergence plots
| File | Purpose |
|------|---------|
| `analyze.R` / `analyze_without_sparsity.r` | General per-parameter error summaries/boxplots (sparse / dense). |
| `convergence_plot.R` | Per-parameter convergence vs `N` — unpenalized, sparse. |
| `convergence_without_sparsity.R` | Per-parameter convergence vs `N` — dense. |
| `convergence_ridge_plot.R` | Per-parameter convergence — fixed-`lambda` ridge, sparse. |
| `convergence_ridge_cv_plot.R` | Per-parameter convergence — CV ridge (random folds), sparse. |
| `convergence_ridge_cv_compare.R` | Compare **random vs site-grouped** CV folds per parameter, sparse. |
| `convergence_ridge_cv_grouped_dense.R` | Per-parameter convergence — site-grouped CV ridge, **dense**. |
| `convergence_theta.R` | Convergence (per-replicate **RMSE**) of the derived potentials `theta_{s,jj}` (main) and `theta_{s,jj'}` (interaction), across all methods and both scenarios. `theta` stays identifiable even when individual params (e.g. `A` in dense) are not. |
| `theta_edge_mean_error.R` | Mean **signed** error of the interaction edge potentials `theta_{s,jj'}` (boxplot + Monte-Carlo line summary). |

### Diagnostics
| File | Purpose |
|------|---------|
| `check_multicollinearity.R` | VIF / correlation structure of the pseudo-likelihood design matrix. |
| `check_sparsity_sparse.R` / `check_sparsity_dense.R` | Response presence-rate / richness diagnostics for each scenario. |
| `check_unpenalized_warnings.R` | Re-run unpenalized fits to capture GLM convergence warnings into `Simulation_Results/warnings/` (does not overwrite estimates). |

## Reproducing

```r
# from the project root, in R
source("Code/FCIR/main.r")                 # sparse scenario  -> Rdata_Sparse/
source("Code/FCIR/main_without_sparsity.r")# dense scenario   -> Rdata_Dense/
source("Code/FCIR/convergence_theta.R")    # figures          -> Simulation_Results/plots/
```

Global grid: `N in {50,100,200,400,800}`, `P in {30,60}`, `L=3`, `K=2`,
`B_reps=1000`, `seed=42`.
