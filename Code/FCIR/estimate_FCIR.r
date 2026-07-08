estimate_unpenalized_FCIR <- function(Y, X, Tr,
                                      standardize = FALSE, returnX_only = FALSE) {
    # Y: N x P binary response matrix
    # X: N x L environment matrix (first column should be 1s for intercept)
    # Tr: P x K species traits matrix
    # standardize: whether to standardize the design matrix (excluding intercept) before fitting the model. All results are then returned on this standardize scale
    # returnX_only: if TRUE, return the design matrix without fitting the model
    
    N = nrow(Y)
    P = ncol(Y)
    L = ncol(X)
    K = ncol(Tr)
    
    # 1. Pre-compute Trait Differences
    Delta = array(0, dim = c(P, P, K))
    for(j in 1:P) {
        for(j_prime in 1:P) {
            Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
        }
    }
    n_obs = N * P               
    n_params = 2*L + 2*L*K     
    
    glm_Y = numeric(n_obs)
    glm_X = matrix(0, nrow = n_obs, ncol = n_params)
    
    row_idx = 1
    for(s in 1:N){
        x_s = X[s, ]     
        y_s = Y[s, ]     
        R_s = sum(y_s)  
        
        for(j in 1:P){
            glm_Y[row_idx] = y_s[j]
            
            # --- Main Effects ---
            t_j = Tr[j, ]
            comp1_beta0 = x_s
            comp2_B     = kronecker(t_j, x_s) # Equivalent to t_j \otimes x_s
            
            # --- Interaction Effects ---
            neighbor_sum = R_s - y_s[j]       # Equivalent to sum_{j' \neq j} y_sj'
            comp3_alpha0 = neighbor_sum * x_s
            
            # Calculate the neighbor trait difference weighted sum (w_sj) for focal species j
            w_sj = numeric(K)
            for(j_prime in 1:P){
                if(j_prime != j && y_s[j_prime] == 1){
                    w_sj = w_sj + Delta[j, j_prime, ]
                }
            }
            comp4_A = kronecker(w_sj, x_s)    # Equivalent to w_sj \otimes x_s
            
            # Concatenate the four blocks into one row for the current observation
            glm_X[row_idx, ] = c(comp1_beta0, comp2_B, comp3_alpha0, comp4_A)
            
            row_idx = row_idx + 1
        }
    }
    
    getsds <- apply(glm_X[,-1], 2, sd)
    if(standardize){
        glm_X[,-1] <- scale(glm_X[,-1])
        }
    
    if(returnX_only) 
        return(glm_X)
    
    
    # 4. Fit Unpenalized Logistic Regression, then extract and reshape parameters
    logistic_reg = glm(y ~ ., data = data.frame(y = glm_Y, glm_X[,-1]), family = binomial)
    # logistic_reg = glm(glm_Y ~ glm_X + 0, family = binomial)
    
    est_coefs = logistic_reg$coefficients
    
    idx = 1
    hat_beta_0  = est_coefs[idx:(idx + L - 1)]; idx = idx + L
    hat_B_vec   = est_coefs[idx:(idx + L*K - 1)]; idx = idx + L*K
    hat_alpha_0 = est_coefs[idx:(idx + L - 1)]; idx = idx + L
    hat_A_vec   = est_coefs[idx:(idx + L*K - 1)]
    
    # Reshape vectors back into L x K fourth-corner matrices
    hat_B_mat = matrix(hat_B_vec, nrow = L, ncol = K)
    hat_A_mat = matrix(hat_A_vec, nrow = L, ncol = K)
    
    return(list(
        beta_0 = hat_beta_0,
        B_mat = hat_B_mat,
        alpha_0 = hat_alpha_0,
        A_mat = hat_A_mat,
        standardize = standardize,
        getsds = getsds,
        glm_model = logistic_reg # Return the original model to allow checking p-values or summary
    ))
    }
