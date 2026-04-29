source("generate_FCIR.r")
source("estimate_FCIR.r")

# 1. Define the parameter grid you want to sweep across
Ns = c(50, 100, 200, 400, 800)
Ps = c(30, 60)

# Global fixed parameters
L = 3          
K = 2          
B_reps = 1000  
seed = 42      

# Record the total start time
total_start = Sys.time()

# 2. Start the nested loops
for (n in Ns) {
  for (p in Ps) {
    
    data_filename = paste0("Simulation_Results/FCIR_data_N", n, "_P", p, ".Rdata")
    est_filename  = paste0("Simulation_Results/FCIR_estimates_N", n, "_P", p, ".Rdata")
    
    if (!dir.exists(dirname(data_filename))) {
      dir.create(dirname(data_filename), recursive = TRUE)
    }
    
    print(paste("Generating data ( N =", n, ", P =", p, ")..."))
    generate_dense_fcir_data(N = n, P = p, L = L, K = K, B_reps = B_reps, 
                             seed = seed, filename = data_filename)
    
    print("Loading data and starting 1000 Monte Carlo estimations...")
    load(data_filename) 
    
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
  }
}
