project_root = "E:/RStudio/thesis"
fcir_code_dir = file.path(project_root, "Code", "FCIR")
rdata_dir = file.path(project_root, "Simulation_Results", "Rdata_Sparse")

source(file.path(fcir_code_dir, "generate_FCIR.r"))
source(file.path(fcir_code_dir, "estimate_FCIR.r"))
source(file.path(fcir_code_dir, "estimate_penalized_FCIR.r"))

# 1. Define the parameter grid you want to sweep across
Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

# Global fixed parameters
L = 3          
K = 2          
B_reps = 1000  
seed = 42      
ridge_lambdas = c(0.5, 1, 2)

format_lambda_label <- function(lambda_value) {
  gsub("\\.", "p", as.character(lambda_value))
}

# Record the total start time
total_start = Sys.time()

# 2. Start the nested loops
for (n in Ns) {
  for (p in Ps) {
    
    data_filename = file.path(rdata_dir, paste0("FCIR_data_N", n, "_P", p, ".Rdata"))
    est_filename  = file.path(rdata_dir, paste0("FCIR_estimates_N", n, "_P", p, ".Rdata"))
    
    if (!dir.exists(dirname(data_filename))) {
      dir.create(dirname(data_filename), recursive = TRUE)
    }
    
    if (!file.exists(data_filename)) {
      print(paste("Generating data ( N =", n, ", P =", p, ")..."))
      generate_dense_fcir_data(N = n, P = p, L = L, K = K, B_reps = B_reps, 
                               seed = seed, filename = data_filename)
    } else {
      print(paste("Data already exists for N =", n, ", P =", p, "- Skipping generation."))
    }
    
    print("Loading data...")
    load(data_filename) 
    
    # ---------------- Unpenalized Estimation ----------------
    if (!file.exists(est_filename)) {
      print(paste("Starting 1000 Monte Carlo Unpenalized estimations for N =", n, "P =", p, "..."))
      est_beta_0  = matrix(NA, nrow = B_reps, ncol = L)
      est_B_mat   = array(NA, dim = c(L, K, B_reps))
      est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
      est_A_mat   = array(NA, dim = c(L, K, B_reps))
      
      for (b in 1:B_reps) {
        Y_b = Y[, , b] 
        result = estimate_unpenalized_FCIR(Y = Y_b, X = X, Tr = Tr)
        
        est_beta_0[b, ]  = result$beta_0
        est_B_mat[,, b]  = result$B_mat
        est_alpha_0[b, ] = result$alpha_0
        est_A_mat[,, b]  = result$A_mat
      }
      
      save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat,
           beta_0, B_mat, alpha_0, A_mat,
           N = n, P = p, L, K, B_reps, seed,
           file = est_filename)
    } else {
      print(paste("Unpenalized estimates already exist for N =", n, ", P =", p, "- Skipping."))
    }
    
    # ---------------- Ridge Estimation with Fixed Lambda ----------------
    for (lambda_value in ridge_lambdas) {
      lambda_label = format_lambda_label(lambda_value)
      est_filename_ridge = file.path(
        rdata_dir,
        paste0("FCIR_estimates_ridge_lambda", lambda_label, "_N", n, "_P", p, ".Rdata")
      )
      
      if (!file.exists(est_filename_ridge)) {
        print(paste("Starting Ridge estimations for N =", n, "P =", p, "lambda =", lambda_value, "..."))
        est_beta_0  = matrix(NA, nrow = B_reps, ncol = L)
        est_B_mat   = array(NA, dim = c(L, K, B_reps))
        est_alpha_0 = matrix(NA, nrow = B_reps, ncol = L)
        est_A_mat   = array(NA, dim = c(L, K, B_reps))
        
        for (b in 1:B_reps) {
          Y_b = Y[, , b] 
          result_ridge = estimate_penalized_FCIR(
            Y = Y_b,
            X = X,
            Tr = Tr,
            alpha = 0,
            lambda = lambda_value,
            use_cv = FALSE
          )
          
          est_beta_0[b, ]  = result_ridge$beta_0
          est_B_mat[,, b]  = result_ridge$B_mat
          est_alpha_0[b, ] = result_ridge$alpha_0
          est_A_mat[,, b]  = result_ridge$A_mat
        }
        
        penalty = "ridge"
        alpha_value = 0
        save(est_beta_0, est_B_mat, est_alpha_0, est_A_mat,
             beta_0, B_mat, alpha_0, A_mat,
             N = n, P = p, L, K, B_reps, seed,
             penalty, alpha_value, lambda_value,
             file = est_filename_ridge)
      } else {
        print(paste("Ridge estimates already exist for N =", n, ", P =", p,
                    ", lambda =", lambda_value, "- Skipping."))
      }
    }
    
  }
}

print(paste("Total execution time:", Sys.time() - total_start))
