library(car)
library(corrplot)

project_root <- "F:/ising model thesis/-Fourth-Corner-Ising-Model-"
rdata_dir <- file.path(project_root, "Simulation_Results", "FCIR", "Rdata_Sparse")
plot_dir <- file.path(project_root, "Simulation_Results", "FCIR", "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

Ns <- c(50, 100, 200, 400, 800)
Ps <- c(30, 60)

for (n in Ns) {
  for (p in Ps) {
    
    # Specify the data file to read
    data_file <- file.path(rdata_dir, paste0("FCIR_data_N", n, "_P", p, ".Rdata"))
    
    if (!file.exists(data_file)) {
      warning(paste("Data file not found:", data_file, "- Skipping."))
      next
    }
    
    load(data_file)
    cat("Successfully loaded data: N =", n, ", P =", p, ", L =", L, ", K =", K, "\n\n")
    
    # Extract the Y generated from the first Monte-Carlo simulation
    Y_b <- Y[, , 1]
    
    # ---------------------------------------------------------
    # 1. Reconstruct the pseudo-likelihood design matrix glm_X
    # (Logic is identical to estimate_unpenalized_FCIR.r)
    # ---------------------------------------------------------
    Delta <- array(0, dim = c(P, P, K))
    for(j in 1:P) {
      for(j_prime in 1:P) {
        Delta[j, j_prime, ] <- abs(Tr[j, ] - Tr[j_prime, ])
      }
    }
    
    n_obs <- n * p
    n_params <- 2*L + 2*L*K
    glm_Y <- numeric(n_obs)
    glm_X <- matrix(0, nrow = n_obs, ncol = n_params)
    
    row_idx <- 1
    for(s in 1:n){
      x_s <- X[s, ]     
      y_s <- Y_b[s, ]     
      R_s <- sum(y_s)   
      
      for(j in 1:p){
        glm_Y[row_idx] <- y_s[j]
        
        t_j <- Tr[j, ]
        comp1_beta0 <- x_s
        comp2_B     <- kronecker(t_j, x_s) 
        
        neighbor_sum <- R_s - y_s[j]       
        comp3_alpha0 <- neighbor_sum * x_s
        
        w_sj <- numeric(K)
        for(j_prime in 1:p){
          if(j_prime != j && y_s[j_prime] == 1){
            w_sj <- w_sj + Delta[j, j_prime, ]
          }
        }
        comp4_A <- kronecker(w_sj, x_s)    
        
        glm_X[row_idx, ] <- c(comp1_beta0, comp2_B, comp3_alpha0, comp4_A)
        row_idx <- row_idx + 1
      }
    }
    

    names_beta0 <- paste0("X_s_", 1:L)
    

    names_B <- paste0("t_j_", rep(1:K, each=L), "_X_s_", rep(1:L, times=K))
    

    names_alpha0 <- paste0("R_s_X_s_", 1:L)
    

    names_A <- paste0("w_sj_", rep(1:K, each=L), "_X_s_", rep(1:L, times=K))
    
    col_names <- c(names_beta0, names_B, names_alpha0, names_A)
    colnames(glm_X) <- col_names
    
    
    # ---------------------------------------------------------
    # 2. Collinearity Detection: Condition Number
    # ---------------------------------------------------------
    # The first column of glm_X is typically the global intercept (since X_s_1 = 1)
    # We exclude this constant intercept term when calculating correlation matrix and VIF
    glm_X_no_intercept <- glm_X[, -1]
    
    # Calculate the standardized correlation matrix of the design matrix
    cor_matrix <- cor(glm_X_no_intercept)
    
    # Calculate eigenvalues and Condition Number (Kappa)
    eigenvalues <- eigen(cor_matrix)$values
    condition_number <- max(eigenvalues) / min(eigenvalues)
    kappa_index <- sqrt(condition_number)
    
    cat("Condition Number:", round(condition_number, 2), "\n")
    cat("Kappa Index (sqrt of Cond Num):", round(kappa_index, 2), "\n")
    
    plot_filename <- file.path(plot_dir, paste0("Correlation_DesignMatrix_glmX_N", n, "_P", p, ".png"))
    png(plot_filename, width = 1000, height = 1000, res = 120)
    corrplot(cor_matrix, method="color", type="upper", tl.col="black", tl.cex=0.8,
             title = paste("Correlation within the Design Matrix (N =", n, ", P =", p, ")"), 
             mar = c(0,0,2,0))
    dev.off()
    
    
  }
}
