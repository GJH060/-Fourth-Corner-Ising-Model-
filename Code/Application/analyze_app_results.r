# Post-processing of the three Southern Ocean CPR fits.
#
#   Rscript Code/Application/analyze_app_results.r
#   source("Code/Application/analyze_app_results.r")
#
# Reads the three parameter blocks saved by main_adaptive_lasso_FCIR{,_I,_M}.r
# and turns them into labelled tables and figures. Nothing here refits anything.
#
# Two structural facts drive how these are laid out, and both come from Traitmat
# being the full indicator expansion of broad_func_grp:
#
#   (1) The trait blocks sum to the beta_0 block, so glm() aliases the last trait
#       column and Pteropod acts as the reference group -- B[, K] is exactly 0 by
#       construction in FCIR and FCIR_M, not by selection. It is dropped from the
#       B heatmaps and flagged in the tables.
#
#   (2) Delta_jj' = |t_j - t_j'| is the zero vector for every same-group pair
#       (verified: exactly the 99 same-group pairs of 435). So the fourth-corner
#       interaction A cannot see within-group association at all -- all 99 of
#       those pairs are absorbed into alpha_0. FCIR_M's free Theta_int can see
#       them, which makes the within- vs cross-group comparison of Theta_int a
#       direct check on whether A's parameterisation is losing real signal.

library(ggplot2)
library(dplyr)
library(tidyr)

project_root = "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
app_data_dir = file.path(project_root, "Application")
rdata_dir = file.path(project_root, "Simulation_Results", "Application", "Rdata")
tab_dir = file.path(project_root, "Simulation_Results", "Application", "tables")
plot_dir = file.path(project_root, "Simulation_Results", "Application", "plots")
for (d in c(tab_dir, plot_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

tol = 1e-8

# --- Data and labels --------------------------------------------------------
load(file.path(app_data_dir, "reduceddat.RData"))   # resp, X, Traitmat

Y = as.matrix(resp); storage.mode(Y) = "numeric"
Tr = Traitmat
P = ncol(Y); K = ncol(Tr)

species = colnames(Y)
group = sub("^broad_func_grp", "", colnames(Tr))[apply(Tr, 1, function(r) which(r == 1))]
grp_levels = sub("^broad_func_grp", "", colnames(Tr))
env_labels = c("(const)", "salinity", "water_temp", "PAR")

prevalence = colMeans(Y)
# Species ordered by functional group, then by prevalence within group -- used
# for every species-indexed figure so rows and columns stay comparable.
sp_order = order(match(group, grp_levels), -prevalence)

f_fcir = file.path(rdata_dir, "app_FCIR_adaptive_lasso_unpenalized_N29983.Rdata")
f_fcir_i = file.path(rdata_dir, "app_FCIR_I_adaptive_lasso_unpenalized_N29983.Rdata")
f_fcir_m = file.path(rdata_dir, "app_FCIR_M_adaptive_lasso_unpenalized_N29983.Rdata")
stopifnot(file.exists(f_fcir), file.exists(f_fcir_i), file.exists(f_fcir_m))

e_f = new.env(); load(f_fcir, envir = e_f)
e_i = new.env(); load(f_fcir_i, envir = e_i)
e_m = new.env(); load(f_fcir_m, envir = e_m)

stopifnot(identical(species, rownames(Tr)),
          e_f$P == P, e_i$P == P, e_m$P == P, e_f$L == length(env_labels))

label_mat <- function(M) { dimnames(M) = list(env_labels, grp_levels); M }
e_f$est_B_mat = label_mat(e_f$est_B_mat); e_f$est_A_mat = label_mat(e_f$est_A_mat)
e_i$est_A_mat = label_mat(e_i$est_A_mat); e_m$est_B_mat = label_mat(e_m$est_B_mat)

# --- Shared plotting helpers ------------------------------------------------
# Diverging fill centred at zero, symmetric limits so +x and -x read equally.
# `clip` squishes values beyond the limit onto the end colour instead of dropping
# them: Theta_int carries two separation artefacts at -7.6 and -4.7 that would
# otherwise flatten every real edge to near-white.
fill_div <- function(vals, name, clip = NULL) {
  lim <- if (is.null(clip)) max(abs(vals), na.rm = TRUE) else clip
  scale_fill_gradient2(name = name, low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1) * lim,
                       oob = scales::squish, na.value = "grey88")
}

theme_hm = theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

# Long form of an L x K fourth-corner matrix, tagged with the model name.
tidy_lk <- function(M, model) {
  as.data.frame(as.table(M), stringsAsFactors = FALSE) |>
    setNames(c("env", "group", "estimate")) |>
    mutate(model = model,
           env = factor(env, levels = env_labels),
           group = factor(group, levels = grp_levels),
           selected = abs(estimate) > tol)
}

# Heatmap of fourth-corner blocks; zeros (screened out by the lasso) are left
# blank rather than painted white, so "shrunk to zero" is visually distinct
# from "estimated near zero".
hm_lk <- function(df, title, subtitle) {
  ggplot(df, aes(group, env)) +
    geom_tile(aes(fill = ifelse(selected, estimate, NA)), colour = "grey85") +
    geom_text(aes(label = ifelse(selected, sprintf("%.2f", estimate), "0")),
              size = 3, colour = "grey20") +
    fill_div(df$estimate, "estimate") +
    facet_wrap(~ model, ncol = 1) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_hm
}

# --- 1. Model overview ------------------------------------------------------
overview = data.frame(
  model = c("FCIR", "FCIR_I", "FCIR_M"),
  main_effect = c("fourth-corner B", "free beta_j (P x L)", "fourth-corner B"),
  interaction = c("fourth-corner A", "fourth-corner A", "free Theta_int (435 edges)"),
  n_params = c(length(e_f$penalty_factor), length(e_i$penalty_factor), length(e_m$penalty_factor)),
  n_nonzero = c(sum(abs(c(e_f$est_beta_0, e_f$est_B_mat, e_f$est_alpha_0, e_f$est_A_mat)) > tol),
                sum(abs(c(e_i$est_Beta_mat, e_i$est_alpha_0, e_i$est_A_mat)) > tol),
                sum(abs(c(e_m$est_beta_0, e_m$est_B_mat,
                          e_m$est_Theta_int[upper.tri(e_m$est_Theta_int)])) > tol)),
  lambda_min = c(e_f$est_lambda, e_i$est_lambda, e_m$est_lambda),
  glm_aliased = c(e_f$init_aliased, e_i$init_aliased, e_m$init_aliased),
  elapsed_s = round(c(e_f$elapsed, e_i$elapsed, e_m$elapsed), 1),
  quasi_separation = c(any(grepl("0 or 1|\u96f6\u6216\u4e00", e_f$fit_warnings)),
                       any(grepl("0 or 1|\u96f6\u6216\u4e00", e_i$fit_warnings)),
                       any(grepl("0 or 1|\u96f6\u6216\u4e00", e_m$fit_warnings))),
  stringsAsFactors = FALSE)
overview$pct_nonzero = round(100 * overview$n_nonzero / overview$n_params, 1)
write.csv(overview, file.path(tab_dir, "app_model_overview.csv"), row.names = FALSE)
print(overview)

# --- 2. Fourth-corner main effect B: FCIR vs FCIR_M -------------------------
# Pteropod is the aliased reference group in both models: its column is zero by
# construction, so it is reported in the table but dropped from the figure.
B_long = bind_rows(tidy_lk(e_f$est_B_mat, "FCIR"), tidy_lk(e_m$est_B_mat, "FCIR_M")) |>
  mutate(is_reference_group = group == grp_levels[K],
         group_size = as.integer(table(group)[as.character(group)]))
write.csv(B_long, file.path(tab_dir, "app_B_fourth_corner.csv"), row.names = FALSE)

ggsave(file.path(plot_dir, "app_B_fourth_corner.png"),
       hm_lk(filter(B_long, !is_reference_group),
             "Fourth-corner main effect B",
             sprintf("environment x functional group; %s is the aliased reference group (B[, K] = 0 by construction)",
                     grp_levels[K])),
       width = 8, height = 6, dpi = 300)

# --- 3. Fourth-corner interaction A: FCIR vs FCIR_I -------------------------
# No aliasing here: A is identified in both models (FCIR_I has no beta_0 block
# at all, and in FCIR the alias falls on B).
A_long = bind_rows(tidy_lk(e_f$est_A_mat, "FCIR"), tidy_lk(e_i$est_A_mat, "FCIR_I"))
write.csv(A_long, file.path(tab_dir, "app_A_fourth_corner.csv"), row.names = FALSE)

ggsave(file.path(plot_dir, "app_A_fourth_corner.png"),
       hm_lk(A_long, "Fourth-corner interaction A",
             "cross-group association only: Delta = 0 for all 99 same-group pairs, which alpha_0 absorbs"),
       width = 8, height = 6, dpi = 300)

# --- 4. FCIR_I species-specific main effects --------------------------------
Beta_long = as.data.frame(as.table(structure(e_i$est_Beta_mat,
                                             dimnames = list(species, env_labels))),
                          stringsAsFactors = FALSE) |>
  setNames(c("species", "env", "estimate")) |>
  mutate(group = group[match(species, species)],
         prevalence = prevalence[species],
         selected = abs(estimate) > tol,
         env = factor(env, levels = env_labels),
         species = factor(species, levels = species[sp_order]))
write.csv(Beta_long, file.path(tab_dir, "app_FCIR_I_Beta.csv"), row.names = FALSE)

# The intercept column and the three environmental slopes differ by an order of
# magnitude (beta_j[, 1] spans -6..1, the slopes -0.4..0.5), so a shared fill
# scale flattens every slope to white. They are plotted separately: the slopes
# as a heatmap, the intercept -- which is one number per species -- as a bar.
# FCIR_I selected lambda.min = 1.2e-4, effectively no penalty, and its glm
# initialiser hit quasi-separation -- so a few slopes run away (metridia_lucens
# on PAR reaches -4.4 against a 98th percentile of ~0.4). Clip here too.
beta_slopes = filter(Beta_long, env != "(const)")
beta_clip = round(as.numeric(quantile(abs(beta_slopes$estimate[beta_slopes$selected]), 0.98)), 2)

g_beta = ggplot(beta_slopes, aes(env, species)) +
  geom_tile(aes(fill = ifelse(selected, estimate, NA)), colour = "grey90") +
  fill_div(beta_slopes$estimate, "beta_j", clip = beta_clip) +
  labs(title = "FCIR_I species-specific environmental responses",
       subtitle = sprintf("slopes only; grey = exactly zero
scale clipped at +/-%.2f (%d of %d exceed it, min %.1f)
species ordered by functional group, then prevalence",
                          beta_clip, sum(abs(beta_slopes$estimate) > beta_clip),
                          nrow(beta_slopes), min(beta_slopes$estimate)),
       x = NULL, y = NULL) + theme_hm
ggsave(file.path(plot_dir, "app_FCIR_I_Beta.png"), g_beta, width = 7.5, height = 8, dpi = 300)

g_beta0 = ggplot(filter(Beta_long, env == "(const)"), aes(estimate, species)) +
  geom_col(aes(fill = group), width = 0.75) +
  geom_vline(xintercept = 0, colour = "grey40") +
  labs(title = "FCIR_I species intercepts beta_j[1]",
       subtitle = "baseline log-odds of presence at mean environment; same ordering as the slope heatmap",
       x = "beta_j[1]", y = NULL, fill = "functional group") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "right")
ggsave(file.path(plot_dir, "app_FCIR_I_Beta_intercepts.png"), g_beta0,
       width = 8, height = 8, dpi = 300)

# --- 5. FCIR_M free association matrix Theta_int ----------------------------
Theta = e_m$est_Theta_int
dimnames(Theta) = list(species, species)
ut = which(upper.tri(Theta), arr.ind = TRUE)

# Co-occurrence support for each pair. An edge estimated from a handful of joint
# presences is the signature of quasi-separation, not of a strong interaction,
# so n11 is what separates a real edge from a numerical artefact.
n11 = crossprod(Y)          # P x P, n11[j, j'] = joint presences
n_j = colSums(Y)

edges = data.frame(
  sp1 = species[ut[, 1]], sp2 = species[ut[, 2]],
  grp1 = group[ut[, 1]], grp2 = group[ut[, 2]],
  theta = Theta[ut],
  prev1 = prevalence[ut[, 1]], prev2 = prevalence[ut[, 2]],
  n11 = n11[ut], n10 = n_j[ut[, 1]] - n11[ut], n01 = n_j[ut[, 2]] - n11[ut],
  stringsAsFactors = FALSE) |>
  mutate(same_group = grp1 == grp2,
         n00 = nrow(Y) - n11 - n10 - n01,
         min_cell = pmin(n11, n10, n01, n00),   # 0 => complete separation of the 2x2 table
         selected = abs(theta) > tol) |>
  arrange(desc(abs(theta)))
write.csv(edges, file.path(tab_dir, "app_FCIR_M_Theta_edges.csv"), row.names = FALSE)

cat("\n--- Theta_int: 15 strongest edges ---\n")
print(head(edges[, c("sp1", "sp2", "same_group", "theta", "n11", "min_cell")], 15))

# Theta_int heatmap, species blocked by functional group. The colour scale is
# clipped at the 98th percentile of |theta| so the bulk of the matrix stays
# readable; the handful above it are the separation artefacts of section 8.
theta_clip = round(as.numeric(quantile(abs(edges$theta[edges$selected]), 0.98)), 1)

# Theta_int heatmap, species blocked by functional group.
Theta_long = as.data.frame(as.table(Theta), stringsAsFactors = FALSE) |>
  setNames(c("sp1", "sp2", "theta")) |>
  filter(sp1 != sp2) |>
  mutate(sp1 = factor(sp1, levels = species[sp_order]),
         sp2 = factor(sp2, levels = species[sp_order]))
grp_breaks = cumsum(table(factor(group[sp_order], levels = grp_levels)))
grp_breaks = grp_breaks[-length(grp_breaks)] + 0.5

g_theta = ggplot(Theta_long, aes(sp1, sp2)) +
  geom_tile(aes(fill = ifelse(abs(theta) > tol, theta, NA))) +
  geom_vline(xintercept = grp_breaks, colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = grp_breaks, colour = "grey40", linewidth = 0.3) +
  fill_div(Theta_long$theta, "theta", clip = theta_clip) +
  labs(title = "FCIR_M residual association matrix Theta_int",
       subtitle = sprintf("%d of %d edges selected; grey = shrunk to exactly zero; lines separate functional groups
colour scale clipped at +/-%.1f; %d edges exceed it, incl. the %d with zero co-occurrence (min theta = %.1f)",
                          sum(edges$selected), nrow(edges), theta_clip,
                          sum(abs(edges$theta) > theta_clip), sum(edges$n11 == 0),
                          min(edges$theta)),
       x = NULL, y = NULL) + theme_hm +
  theme(axis.text = element_text(size = 6))
ggsave(file.path(plot_dir, "app_FCIR_M_Theta_heatmap.png"), g_theta,
       width = 9, height = 8, dpi = 300)

# --- 6. Diagnostic: are the extreme edges quasi-separation artefacts? -------
g_diag = ggplot(filter(edges, selected), aes(n11 + 1, theta)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(aes(colour = min_cell == 0), alpha = 0.7, size = 1.6) +
  scale_x_log10() +
  scale_colour_manual(name = "2x2 table has an empty cell",
                      values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
  labs(title = "Theta_int edge strength against co-occurrence support",
       subtitle = "edges resting on few joint presences are where quasi-separation inflates |theta|",
       x = "joint presences n11 (log scale, +1)", y = "theta") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "app_FCIR_M_edge_support.png"), g_diag,
       width = 7, height = 5, dpi = 300)

# --- 7. Diagnostic: does A's same-group blind spot lose signal? -------------
# FCIR's A is forced to zero on same-group pairs (Delta = 0). FCIR_M is not.
# If same-group edges in Theta_int are systematically as strong as cross-group
# ones, then the fourth-corner interaction is discarding real structure.
sg = edges |> group_by(same_group) |>
  summarise(n_pairs = n(), n_selected = sum(selected),
            pct_selected = round(100 * mean(selected), 1),
            median_abs_theta = round(median(abs(theta)), 3),
            mean_theta = round(mean(theta), 3), .groups = "drop")
write.csv(sg, file.path(tab_dir, "app_Theta_same_vs_cross_group.csv"), row.names = FALSE)
cat("\n--- Theta_int: within- vs cross-functional-group edges ---\n"); print(as.data.frame(sg))

g_sg = ggplot(filter(edges, selected),
              aes(ifelse(same_group, "within group", "cross group"), theta)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_boxplot(outlier.size = 0.8, width = 0.5, fill = "grey92") +
  labs(title = "Theta_int by functional-group membership",
       subtitle = "FCIR's A is structurally zero on within-group pairs;
FCIR_M is free to fit them",
       x = NULL, y = "theta") +
  theme_minimal(base_size = 11)
ggsave(file.path(plot_dir, "app_Theta_same_vs_cross_group.png"), g_sg,
       width = 6, height = 5, dpi = 300)

# --- 8. Diagnostic: recording-protocol turnover masquerading as association --
# The two strongest negative edges (theta = -7.6, -4.7) both have n11 = 0. They
# are not ecology. ctenocalanus_sp_cf_citer is recorded 1991-2014 and
# ctenocalanus_sp_cf_citer_2 only 2011-2016 -- the same taxon under two codes
# either side of a taxonomic revision, so they can never co-occur. Several other
# taxa are likewise confined to a subset of the 24 survey years, and any two
# species whose recording windows do not overlap are pushed toward a strong
# negative theta by construction.
#
# This is protocol heterogeneity, not the spatio-temporal correlation that is
# deliberately out of scope: no autocorrelation model would remove it, because
# those zeros are absences of *recording*, not absences of the organism.
yr = X$year
yr_present = lapply(species, function(s) unique(yr[Y[, s] == 1]))
names(yr_present) = species
n_years_total = length(unique(yr))

edges$years_shared = mapply(function(a, b) length(intersect(yr_present[[a]], yr_present[[b]])),
                            edges$sp1, edges$sp2)
edges$years_jaccard = round(edges$years_shared /
  mapply(function(a, b) length(union(yr_present[[a]], yr_present[[b]])),
         edges$sp1, edges$sp2), 3)
write.csv(edges, file.path(tab_dir, "app_FCIR_M_Theta_edges.csv"), row.names = FALSE)

yr_tab = data.frame(species = species, group = group,
                    prevalence = round(prevalence, 4),
                    n_years_recorded = lengths(yr_present)[species],
                    first_year = sapply(yr_present, min)[species],
                    last_year = sapply(yr_present, max)[species],
                    row.names = NULL) |> arrange(n_years_recorded)
write.csv(yr_tab, file.path(tab_dir, "app_species_recording_windows.csv"), row.names = FALSE)
cat("\n--- species recorded in fewest of the", n_years_total, "survey years ---\n")
print(head(yr_tab, 8))

g_yr = ggplot(filter(edges, selected), aes(years_jaccard, theta)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(aes(colour = n11 == 0), alpha = 0.7, size = 1.6) +
  scale_colour_manual(name = "never co-occurs", values = c("FALSE" = "grey35", "TRUE" = "#B2182B")) +
  labs(title = "Theta_int against overlap of the two species' recording windows",
       subtitle = "pairs whose survey years barely overlap are driven negative by protocol turnover, not ecology",
       x = "Jaccard overlap of years in which each species was recorded", y = "theta") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(plot_dir, "app_Theta_recording_overlap.png"), g_yr,
       width = 7, height = 5, dpi = 300)

cat("\nSpearman correlation(theta, year-window overlap) among selected edges: ",
    round(cor(edges$theta[edges$selected], edges$years_jaccard[edges$selected],
              method = "spearman"), 3), "\n", sep = "")

cat("\nTables ->", tab_dir, "\nPlots  ->", plot_dir, "\n")
