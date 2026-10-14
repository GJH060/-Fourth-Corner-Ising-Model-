# Fourth-Corner Ising Regression Model

Code accompanying the thesis *A fourth-corner Ising regression model*. The
fourth-corner Ising regression (FCIR) model is a site-varying Ising model whose
main effects `theta_{s,jj}` and pairwise interactions `theta_{s,jj'}` are driven
by environmental covariates `x_s` and species traits `t_j`. The thesis also
studies two special cases:

- **FCIR_I**: traits enter only through the pairwise interactions.
- **FCIR_M**: traits enter only through the main effects.

All three models are fitted by adaptive lasso penalised maximum pseudo-likelihood,
written as a single stacked logistic regression and fitted with `glmnet`.

## Requirements

- R (>= 4.6), with the packages `glmnet`, `foreach`, `doParallel`,
  `IsingSampler`, `Matrix` and `dplyr`. The estimator comparison figures
  also need `ggplot2` and `patchwork`.
- Windows PowerShell, for the `*.ps1` scripts that draw the Chapter 4
  selection figures (they use `System.Drawing`).

```r
install.packages(c("glmnet", "foreach", "doParallel", "IsingSampler",
                   "Matrix", "dplyr", "ggplot2", "patchwork"))
```

**Paths.** Most scripts set `project_root` to an absolute path on their first
line; change it to the location of your clone before running. The FCIR_M
scripts read the environment variable `FCIR_PROJECT_ROOT` instead.

## Data

- **Simulated data** are not stored in the repository (`*.Rdata` is
  git-ignored). The `main_generate_data*` scripts below regenerate them
  exactly, since every configuration uses `seed = 2026`.
- **SOCPR data** (Chapter 5) come from the Southern Ocean Continuous Plankton
  Recorder Zooplankton Records, Australian Antarctic Data Centre,
  <https://doi.org/10.26179/5ee84f77cc4ec>. The processed subset used in the
  thesis is `Application/reduceddat.RData`, which contains the response matrix,
  the environmental covariates and the functional-group indicator matrix.

**Reproducibility.** Cross-validation folds are assigned at random. The FCIR_M
script sets a seed for each replicate, but the FCIR and FCIR_I scripts do not,
so rerunning those gives Monte Carlo averages that differ slightly from the
thesis.

## Repository layout

```
Application/            processed SOCPR data (input to Chapter 5; git-ignored)
Application_SOCPR/
  plots/                Chapter 5 figures
  tables/               Chapter 5 tables (CSV, and ASCII .txt)
Code/
  FCIR/                 full FCIR model: data generation, estimation, evaluation
  FCIR_I/               FCIR_I special case
  FCIR_M/               FCIR_M special case
  Application/          Chapter 5 model fits, figures and tables
  *.r                   estimator comparisons shared across the three models
Simulation_Results/
  FCIR/, FCIR_I/, FCIR_M/
    Rdata*/             simulated data and estimates (git-ignored)
    tables/             summary tables (CSV, and ASCII .txt)
    plots/              figures
  Application/Rdata/    Chapter 5 model fits (git-ignored)
  comparison_tables.tex LaTeX bodies of the estimator comparison tables
```

## Thesis figures and the scripts that produce them

Each pipeline runs left to right: generate data, fit, summarise into a table,
then draw the figure. All paths are relative to the repository root.

### Chapter 4: simulation study

| Thesis figure | Figure file (`Simulation_Results/...`) | Plotting script | Input table | Upstream scripts |
|---|---|---|---|---|
| Fig. 4.1, FCIR selection and RMSE | `FCIR/plots/FCIR_sparse_selection_metrics.png`, `FCIR/plots/FCIR_dense_RMSE.png` | `Code/FCIR/plot_FCIR_simulation_results.ps1` | `FCIR/tables/adaptive_lasso_selection_wide_unpenalized_Sparse.csv`, `FCIR/tables/adaptive_lasso_dense_diagnostics_wide_unpenalized_Sparse.csv` | `Code/FCIR/main_generate_data.r` → `Code/FCIR/main_adaptive_lasso.r` → `Code/FCIR/evaluate_selection_adaptive_lasso.r` |
| Fig. 4.2, adaptive vs ordinary lasso selection | `FCIR/plots/FCIR_lasso_vs_adaptive_selection.png` | `Code/plot_estimator_comparison.r` | `FCIR/tables/lasso_selection_wide.csv` | `Code/FCIR/main_lasso_and_unpenalized.r` → `Code/compare_estimators_RMSE.r` |
| Fig. 4.3, FCIR RMSE of three estimators | `FCIR/plots/FCIR_RMSE_estimator_comparison.png` | `Code/plot_estimator_comparison.r` | `FCIR/tables/estimator_RMSE_comparison.csv` | `Code/FCIR/main_lasso_and_unpenalized.r` → `Code/compare_estimators_RMSE.r` |
| Fig. 4.4, FCIR_I selection and RMSE | `FCIR_I/plots/FCIR_I_sparse_selection_metrics.png` | `Code/FCIR_I/plot_FCIR_I_simulation_results.ps1` | `FCIR_I/tables/Rescale_FCIR_I_adaptive_lasso_selection_wide_sparse_Beta_unpenalized.csv` | `Code/FCIR_I/main_generate_data_FCIR_I.r` (`use_sparse_beta = TRUE`) → `Code/FCIR_I/main_adaptive_lasso_sparse_Beta_FCIR_I.r` → `Code/FCIR_I/evaluate_selection_adaptive_lasso_sparse_Beta_FCIR_I.r` |
| Fig. 4.5, FCIR_I adaptive lasso vs unpenalised | `FCIR_I/plots/FCIR_I_RMSE_estimator_comparison.png` | `Code/plot_estimator_comparison.r` | `FCIR_I/tables/estimator_RMSE_comparison.csv` | `Code/FCIR_I/main_unpenalized_FCIR_I.r` → `Code/compare_estimators_RMSE.r` |
| Fig. 4.6, FCIR_M selection for `Theta` and RMSE | `FCIR_M/plots/FCIR_M_Theta_selection_metrics.png`, `FCIR_M/plots/FCIR_M_main_effect_RMSE.png` | `Code/FCIR_M/plot_FCIR_M_simulation_results.ps1` | `FCIR_M/tables/New_FCIR_M_adaptive_lasso_selection_wide_unpenalized.csv` | `Code/FCIR_M/main_generate_data_FCIR_M_L3.r` → `Code/FCIR_M/main_adaptive_lasso_FCIR_M.r` → `Code/FCIR_M/evaluate_selection_adaptive_lasso_FCIR_M.r` |
| Fig. 4.7, FCIR_M adaptive lasso vs unpenalised | `FCIR_M/plots/FCIR_M_RMSE_estimator_comparison.png` | `Code/plot_estimator_comparison.r` | `FCIR_M/tables/estimator_RMSE_comparison.csv` | `Code/FCIR_M/main_unpenalized_FCIR_M_L3.r` → `Code/compare_estimators_RMSE.r` |

Simulation settings: `N` in {50, 100, 200, 400, 800}, `L = 3`, `K = 2`,
1000 replicates per configuration, `seed = 2026`. `P` is {30, 60} for FCIR and
FCIR_I, and {10, 20} for FCIR_M.

**Note on FCIR_M.** The thesis uses `L = 3` with all non-intercept coefficients,
including `B`, penalised. `Code/FCIR_M/main_adaptive_lasso_FCIR_M.r` and
`evaluate_selection_adaptive_lasso_FCIR_M.r` are set to this configuration
(`L = 3`, `unpenalize_B = FALSE`). They write files tagged `L3`, for example
`New_FCIR_M_adaptive_lasso_selection_wide_L3_unpenalized.csv`. The untagged
table listed above holds the results reported in the thesis. A rerun of the
`L3` pipeline reproduces them up to Monte Carlo error in the cross-validation
folds. The `L = 1` / `unpenalize_B = TRUE` setting is a separate experiment
that is not reported in the thesis.

### Chapter 5: application to the SOCPR data

| Thesis figure | Figure file (`Application_SOCPR/plots/...`) |
|---|---|
| Fig. 5.1, `B` under FCIR and FCIR_M | `app_B_fourth_corner.png` |
| Fig. 5.2, `A` under FCIR and FCIR_I | `app_A_fourth_corner.png` |
| Fig. 5.3, FCIR_I species-specific intercepts | `app_FCIR_I_Beta_intercepts.png` |
| Fig. 5.4, FCIR_I species-specific slopes | `app_FCIR_I_Beta.png` |
| Fig. 5.5, FCIR_M `Theta` heatmap | `app_FCIR_M_Theta_heatmap.png` |

Pipeline:

1. `Code/Application/main_adaptive_lasso_FCIR.r`,
   `main_adaptive_lasso_FCIR_I.r` and `main_adaptive_lasso_FCIR_M.r` fit the three
   models to `Application/reduceddat.RData`.
2. They save the fits to `Simulation_Results/Application/Rdata/`.
3. `Code/Application/analyze_app_results.r` produces all Chapter 5 figures and
   tables, and writes them to `Application_SOCPR/`.
4. Optionally, `Code/Application/format_app_tables.r` writes ASCII versions of
   the tables.

## Thesis tables and their source files

| Thesis table | Source file (`Simulation_Results/...` unless stated) |
|---|---|
| Table 5.1, model overview | `Application_SOCPR/tables/app_model_overview.csv` |
| Table 5.2, `beta_0` and `alpha_0` estimates | `Application_SOCPR/tables/app_beta0_alpha0.csv` |
| Table 5.3, `Theta` within vs across groups | `Application_SOCPR/tables/app_Theta_same_vs_cross_group.csv` |
| Appendix: FCIR adaptive lasso selection (`B`, `A`) | `FCIR/tables/adaptive_lasso_selection_wide_unpenalized_Sparse.csv` |
| Appendix: FCIR RMSE of `beta_0`, `alpha_0` | `FCIR/tables/adaptive_lasso_dense_diagnostics_wide_unpenalized_Sparse.csv` |
| Appendix: FCIR ordinary lasso selection | `FCIR/tables/lasso_selection_wide.csv` |
| Appendix: FCIR_I selection | `FCIR_I/tables/Rescale_FCIR_I_adaptive_lasso_selection_wide_sparse_Beta_unpenalized.csv` |
| Appendix: FCIR_M selection and RMSE | `FCIR_M/tables/New_FCIR_M_adaptive_lasso_selection_wide_unpenalized.csv` |
| Appendix: estimator RMSE comparisons (all three models) | `{FCIR,FCIR_I,FCIR_M}/tables/estimator_RMSE_comparison.csv`, typeset into `comparison_tables.tex` by `Code/format_comparison_tables.r` |

## Core functions

| File | Purpose |
|---|---|
| `Code/FCIR/generate_FCIR.r`, `Code/FCIR_I/generate_sparse_Beta_FCIR_I.r`, `Code/FCIR_M/generate_FCIR_M.r` | Simulate presence-absence data from each model with `IsingSampler` |
| `Code/FCIR/estimate_FCIR.r`, `Code/FCIR_I/estimate_FCIR_I.r`, `Code/FCIR_M/estimate_FCIR_M.r` | Build the stacked pseudo-likelihood model matrix and fit the unpenalised estimator |
| `Code/FCIR/estimate_adaptive_lasso_FCIR.r`, `Code/FCIR_I/estimate_adaptive_lasso_FCIR_I.r`, `Code/FCIR_M/estimate_adaptive_lasso_FCIR_M.r` | Standardise the model matrix, build the adaptive weights, choose `lambda` by site-grouped cross-validation (`cv.glmnet`) and back-transform the estimates |

## Other scripts

The following are exploratory and not used for any result in the thesis:

- In `Code/FCIR/`: `main.r`, `main_without_sparsity.r`, the `convergence_*`
  scripts, the `*_fixedN*` scripts and the `*_weaken_alpha0_and_Amat*` scripts.
- The folders `Code/Ising_joint/` and `Code/Ising_symmetrizingapproach_NOTUSED/`.
- The `toyexample*.R` scripts at the repository root.
