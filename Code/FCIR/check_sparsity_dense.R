library(ggplot2)
library(dplyr)
library(tidyr)

# Set combinations
Ns <- c(50, 100, 200, 400, 800)
Ps <- c(30, 60)

project_root <- "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
data_dir <- file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Dense")
out_dir <- file.path(project_root, "Simulation_Results", "FCIR", "plots", "sparsity")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

results_list <- list()
site_results_list <- list()

for (n in Ns) {
  for (p in Ps) {
    file_path <- file.path(data_dir, paste0("FCIR_data_Dense_N", n, "_P", p, ".Rdata"))
    if (file.exists(file_path)) {
      cat("Processing N =", n, ", P =", p, "\n")
      # Load into a new environment to prevent variable conflicts
      env <- new.env()
      load(file_path, envir = env)
      Y <- env$Y  # Dimension is expected to be [N, P, B_reps]
      
      # 1. Overall sparsity across all 1000 repetitions
      overall <- sum(Y == 0) / length(Y)
      
      # 2. Species sparsity: average sparsity for each species across all N sites and 1000 repetitions
      # apply over the 2nd dimension (species)
      species_sparsity <- apply(Y, 2, function(x) sum(x == 0) / length(x))
      
      temp_df <- data.frame(
        N = as.factor(n),
        P = as.factor(p),
        Species = 1:p,
        Sparsity = species_sparsity,
        Overall = overall
      )
      
      results_list[[paste0(n, "_", p)]] <- temp_df
      
      # 3. Site sparsity: average sparsity for each site across all P species and 1000 repetitions
      # apply over the 1st dimension (sites)
      site_sparsity <- apply(Y, 1, function(x) sum(x == 0) / length(x))
      
      temp_site_df <- data.frame(
        N = as.factor(n),
        P = as.factor(p),
        Site = 1:n,
        Sparsity = site_sparsity
      )
      
      site_results_list[[paste0(n, "_", p)]] <- temp_site_df
      
    } else {
      cat("File not found:", file_path, "\n")
    }
  }
}

all_results <- bind_rows(results_list)
all_site_results <- bind_rows(site_results_list)

# Plot 1: Boxplot of Species Sparsity across N and P
p1 <- ggplot(all_results, aes(x = N, y = Sparsity, fill = P)) +
  geom_boxplot() +
  labs(title = "Species Sparsity Distribution (Dense Data)",
       x = "Sample Size (N)",
       y = "Sparsity (Proportion of 0s per Species)",
       fill = "Species Size (P)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(limits = c(0, 1))

plot1_path <- file.path(out_dir, "Sparsity_Species_Boxplot_Dense.png")
ggsave(plot1_path, p1, width = 8, height = 6)
cat("Saved Plot 1:", plot1_path, "\n")

# Plot 2: Overall Sparsity Trend
overall_df <- all_results %>% 
  group_by(N, P) %>% 
  summarise(Overall_Sparsity = first(Overall), .groups = 'drop')

p2 <- ggplot(overall_df, aes(x = as.numeric(as.character(N)), y = Overall_Sparsity, color = P, group = P)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(title = "Overall Sparsity Trend (Dense Data)",
       x = "Sample Size (N)",
       y = "Overall Sparsity",
       color = "Species Size (P)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(limits = c(0, 1))

plot2_path <- file.path(out_dir, "Sparsity_Overall_Trend_Dense.png")
ggsave(plot2_path, p2, width = 8, height = 6)
cat("Saved Plot 2:", plot2_path, "\n")

# Plot 3: Boxplot of Site Sparsity across N and P
p3 <- ggplot(all_site_results, aes(x = N, y = Sparsity, fill = P)) +
  geom_boxplot() +
  labs(title = "Site Sparsity Distribution (Dense Data)",
       x = "Sample Size (N)",
       y = "Sparsity (Proportion of 0s per Site)",
       fill = "Species Size (P)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(limits = c(0, 1))

plot3_path <- file.path(out_dir, "Sparsity_Site_Boxplot_Dense.png")
ggsave(plot3_path, p3, width = 8, height = 6)
cat("Saved Plot 3:", plot3_path, "\n")

cat("Sparsity analysis complete!\n")