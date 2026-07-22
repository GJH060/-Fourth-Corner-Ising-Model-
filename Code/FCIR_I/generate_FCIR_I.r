library(IsingSampler)

generate_fcir_I_data <- function(N, P, L, K, B_reps, seed, filename){
  # Restricted Model: FCIR_I 
  # Traits drive the species interaction network only. 
  # Main effects are unconstrained species-specific environmental responses.
  
    set.seed(seed)
    # Beta_mat: P x L matrix. Each row is beta_j for species j.
    # Parameter space expands from L(K+1) to P*L.
    Beta_mat = matrix(runif(P * L, -1, 1), nrow = P, ncol = L)
    
    # Interaction Effect Parameters.
    alpha_0 = runif(L, -1, 1)*0.4
    A_mat = matrix(runif(L * K, -1, 1)*0.4, nrow = L, ncol = K) 
    A_mat[order(abs(A_mat))[1:3]] = 0
    # Normalize the interaction parameters by (P - 1) so the total neighbor field
    # stays O(1) across P (avoids saturation). The conditional-logit formula below
    # is left unchanged; the normalization lives entirely in the parameters.
    #alpha_0 = alpha_0 / (P - 1)
    #A_mat = A_mat / (P - 1)
    
    
    Y = array(data = NA, dim = c(N, P, B_reps))
    X = matrix(rnorm(N * L), nrow = N, ncol = L)
    X[,1] = 1 
    Tr = matrix(runif(P * K, min = -1, max = 1), nrow = P, ncol = K)
    
    Delta = array(0, dim = c(P, P, K))
    for(j in 1:P) {
        for(j_prime in 1:P) {
            Delta[j, j_prime, ] = abs(Tr[j, ] - Tr[j_prime, ])
        }
    }
    
    Beta_temp = 1 
    
    # Pre-compute the site-specific fields once per site (they do not depend on
    # the replicate index b), matching generate_FCIR.r. This avoids recomputing
    # the O(P^2) interaction network B_reps times per site.
    #   theta_jj_all[s, ] = x_s^T beta_j                       (main effects)
    #   Theta_all[[s]][j, jp] = x_s^T alpha_0 + sum_k (x_s^T A_mat)[k] Delta[j, jp, k]
    theta_jj_all = matrix(0, nrow = N, ncol = P)
    Theta_all = vector("list", N)
    for (s in 1:N) {
        x_s = X[s, ]
        theta_jj_all[s, ] = as.numeric(Beta_mat %*% x_s)
        
        v_s = as.numeric(x_s %*% A_mat)              # length K: (x_s^T A_mat)
        Theta_s = matrix(sum(x_s * alpha_0), nrow = P, ncol = P)
        for (k in 1:K) Theta_s = Theta_s + Delta[, , k] * v_s[k]
        diag(Theta_s) = 0                            # no self-interaction
        Theta_all[[s]] = Theta_s
    }
    
    for(b in 1:B_reps){
        Y_b = matrix(NA, nrow = N, ncol = P)
        for(s in 1:N){
            Y_b[s, ] = IsingSampler(1, Theta_all[[s]], theta_jj_all[s, ], Beta_temp, 1000/P, 
                                    responses = c(0L, 1L), method = "MH")
            }
        Y[,,b] = Y_b
        }
    save(Y, X, Tr, Beta_mat, alpha_0, A_mat, N, P, L, K, B_reps, seed, file = filename)
    }