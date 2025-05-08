###############################################################################################
## NAM CAVI FUNCTION
###############################################################################################
if (!require(Rcpp)) install.packages("Rcpp", dependencies = TRUE); suppressPackageStartupMessages(library(Rcpp))
if (!require(mclust)) install.packages("mclust", dependencies = TRUE); suppressPackageStartupMessages(library(mclust))

sourceCpp("NAM_VI.cpp")
sourceCpp("ELBO.cpp")

Run_CAM_VI <- function(X, Y, H = H, L = H, maxIter = maxIter, epsilon = epsilon,  warmstart = TRUE){
  J = ncol(X) # Number of groups
  q = nrow(X) # Dimension of group-level variables
  p = nrow(Y[[1]]) # Dimension of observation-level variables
  n = unlist(lapply(Y, ncol)) # Sample sizes
  
  m0_y = matrix(0, nrow = p); lambda0_y = 0.01; nu0_y = p + 5; W0_y = diag(p); iW0_y = solve(W0_y)
  m0_x = matrix(0, nrow = q); lambda0_x = 0.01; nu0_x = q + 5; W0_x = diag(q); iW0_x = solve(W0_x)
  
  lambdak_x = matrix(stats::rgamma(H,1,10), ncol = 1)
  nuk_x = matrix(c(q + stats::rgamma(H,1,1)), ncol = 1)
  Wk_x = stats::rWishart(n = H, df = nu0_x, Sigma = W0_x)
  
  if(warmstart){
    # 
    # library(mclust)
    out = Mclust(t(X), verbose = FALSE)
    mk_x = cbind(out$parameters$mean,
                 sapply(1:(H-ncol(out$parameters$mean)), function(j){mvrnorm(n = 1, mu = m0_x, Sigma = (1/lambdak_x[j, ])* Wk_x[,,j])}))
    mk_x = mk_x[ , sample(1:ncol(mk_x))]
  }else{
    mk_x =  sapply(1:H, function(j){MASS::mvrnorm(n = 1, mu = m0_x, Sigma = (1/lambdak_x[j, ])* Wk_x[,,j])}) # This is for random mean start with the draw from prior distribution
  }
  
  
  lambdal_y = matrix(stats::rgamma(L,1,10), ncol = 1)
  nul_y = matrix(c(p + stats::rgamma(L,1,1)), ncol = 1)
  Wl_y = stats::rWishart(n = L, df = nu0_y, Sigma = W0_y)
  ml_y =  sapply(1:L, function(j){MASS::mvrnorm(n = 1, mu = m0_y, Sigma = (1/lambdal_y[j, ])* Wl_y[,,j])}) # This is for random mean start with the draw from prior distribution
  
  s1 = 0.1; s2 = 0.1; b_h_bar = s1/s2
  r1 = 0.1; r2 = 0.1; b_lk_bar = r1/r2
  
  ###############################################################################################
  XI_jil = list()
  for(j in 1:J){
    log.XI_il  <- array(stats::rbeta( n[j] * L, 1, 1),dim = c( n[j], L))
    Z           <- apply(log.XI_il, c(1), function(x) matrixStats::logSumExp(x))
    XI_jil[[j]] <- exp(sapply(1:L, function(qq) log.XI_il[,qq]-Z,simplify = "array"))
  }
  ###############################################################################################
  if(H<J){
    log.RHO_jk <- matrix(stats::rbeta(J*H,1,1),J,H)
    Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
    RHO_jk     <- exp(sapply(1:H, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
  }else if(H==J){
    RHO_jk     <- diag(J)
  }else{
    RHO_jk     <- cbind(diag(J), matrix(0,J,H-J))
  }
  
  elbo = 0
  diff_ELBO = 1
  # time at the beginning
  Time.start = Sys.time()
  
  for(iter in 1:maxIter){
    # iter = 1
    if(iter == 1){
      cat(paste0("Iteration: ", iter, "\n"))
    }
    if(iter %% 100 == 0) {
      cat(paste0("Iteration: ", iter, "\n"))
    }
    
    ###############################################################################################
    # First Update Ulk ~ Beta(a_bar_Ulk, b_bar_Ulk), for k = 1,.., H-1; l = 1,.., L
    ###############################################################################################
    UU = Update_Ulk_cpp(XI_jil = XI_jil,
                        RHO_jk = RHO_jk,
                        b_bar = b_lk_bar,
                        L = L,
                        J = J,
                        H = H)
    
    a_bar_Ulk = UU[ , ,1] 
    b_bar_Ulk = UU[ , ,2]
    ElnOM_lk  = UU[ , ,3]
    
    ###############################################################################################
    # Second Update M_jil ~ Multinomial(XI_jil), for j=1,...,J; i = 1,..,n_j; l = 1,.., L
    ###############################################################################################
    XI_jil = Update_XIjil_cpp(Y_grouped = lapply(Y, t),
                              RHO_jk = RHO_jk,
                              ElnOM_lk = ElnOM_lk,
                              Nj = matrix(n, ncol = 1),
                              ml = as.matrix(ml_y),
                              tl = lambdal_y,
                              cl = nul_y,
                              Dl = Wl_y,
                              
                              L = L,
                              J = J,
                              K = H)
    
    ###############################################################################################
    # Third Update Theta_l^Y  l = 1,.., L
    ###############################################################################################
    THETAl_Y  = Update_THETAl_Y_cpp_mvt(Y_grouped = lapply(Y, t),
                                        XI_jil = XI_jil,
                                        m0 = m0_y,
                                        lambda0 = lambda0_y,
                                        nu0 = nu0_y,
                                        W0 = W0_y,
                                        iW0 = iW0_y)
    
    ml_y = THETAl_Y["ml"]; ml_y = ml_y$ml
    lambdal_y = THETAl_Y["lambdal"]; lambdal_y = lambdal_y$lambdal
    nul_y = THETAl_Y["nul"]; nul_y = nul_y = nul_y$nul
    Wl_y =  THETAl_Y["Wl"]; Wl_y = Wl_y$Wl
    Sl_y =  THETAl_Y["Sl_over_Nl"]; Sl_y = Sl_y$Sl_over_Nl
    Ybar =  THETAl_Y["Ybar"]; Ybar = Ybar$Ybar
    Nl_y =  THETAl_Y["Nl"]; Nl_y = Nl_y$Nl
    
    ###############################################################################################
    # Fourth Update v_k ~ Beta(a_vk, b_k), for k = 1,.., H-1
    ###############################################################################################
    var_par_v = Update_Vk_DP_conc_cpp(b_bar = b_h_bar, RHO_jk = RHO_jk);
    
    a_vk     = var_par_v[, 1];
    b_vk     = var_par_v[, 2];
    E_ln_PIk = var_par_v[, 3];
    
    ###############################################################################################
    # Fifth Update S_j ~ Multinomial(RHO_jk), for j=1,...,J; k = 1,.., H
    ###############################################################################################
    RHO_jk = Update_RHOjk_CAM_cpp(XI_jil = XI_jil,
                                  ElnPI_k = E_ln_PIk,
                                  ElnOM_lk = ElnOM_lk,
                                  L = L,
                                  J = J,
                                  H = H)
    
    ###############################################################################################
    # Sixth Update alpha ~ Gamma(s1, s2)
    ###############################################################################################
    Update_alpha = Update_alpha_concentration_par(a_tilde_Vk = a_vk,
                                                  b_tilde_Vk = b_vk,
                                                  conc_hyper = matrix(0.1, nrow = 2))
    
    b_h_bar = Update_alpha[1,]/Update_alpha[2, ]
    
    ###############################################################################################
    # Seventh Update beta ~ Gamma(r1, r2)
    ###############################################################################################
    Update_beta = Update_beta_concentration_par(a_bar_Ulk = a_bar_Ulk,
                                                b_bar_Ulk = b_bar_Ulk,
                                                conc_hyper = matrix(0.1, nrow = 2),
                                                L = L,
                                                H = H)
    
    b_lk_bar = Update_beta[1,]/Update_beta[2, ]
    
    
    
    ###############################################################################################
    # CALCULATION OF ELBO
    ###############################################################################################
    elbo1 = elbo_diff_v(a_tilde_k = a_vk,
                        b_tilde_k = b_vk,
                        S_concDP = Update_alpha)
    
    
    elbo2 = elbo_diff_u(a_bar_Ulk = a_bar_Ulk,
                        b_bar_Ulk = b_bar_Ulk,
                        R_concDP = Update_beta)
    
    elbo3 = elbo_diff_alpha(conc_hyper = matrix(0.1, nrow = 2),
                            S_concDP = Update_alpha)
    
    elbo4 = elbo_diff_beta(conc_hyper = matrix(0.1, nrow = 2),
                           R_concDP = Update_beta)
    
    elbo5 = elbo_diff_S(RHO_jk = RHO_jk,
                        ElnPI = E_ln_PIk)
    
    
    elbo6 = elbo_diff_M(XI_jil = XI_jil,
                        RHO_jk = RHO_jk,
                        ElnOM_lk = ElnOM_lk)
    
    
    
    elbo7 = elbo_p_Y(XI_ijl = XI_jil,
                     ml = ml_y,
                     tl = lambdal_y,
                     cl = nul_y,
                     Dl = Wl_y,
                     Sl = Sl_y,    
                     Ybar = Ybar,   
                     Nl = Nl_y)
    
    lCpl0_y  = log_Const_prod_gamma(D = p, nu = lambda0_y)
    LlogB0_y = L * logB_wish(W = W0_y , nu = lambda0_y, lCpl0_y);
    
    elbo8 = elbo_p_THETA_Y(m0 = m0_y, t0 = lambda0_y, c0 = nu0_y,
                           D0 = W0_y, iD0 = iW0_y,
                           lCpl0 = lCpl0_y, LlogB0 = LlogB0_y,
                           ml = ml_y,
                           tl = lambdal_y,
                           cl = nul_y, Dl = Wl_y)
    
    
    elbo9 = elbo_q_THETA_Y(ml = ml_y, tl = lambdal_y,
                           cl = nul_y, Dl = Wl_y)
    
    
    elbo_value = elbo1 +
      elbo2 + 
      elbo3 + 
      elbo4 + 
      elbo5 + 
      elbo6 + 
      elbo7 +
      elbo8 + 
      elbo9
    
    ## STORE THE ELBO VALUE
    elbo[iter] = elbo_value
    
    if(iter > 1){
      diff_ELBO = (elbo[iter]-elbo[iter - 1])}
    
    if( (abs(diff_ELBO) > epsilon) & (iter == maxIter) ) {
      cat(paste0("Warning! Maximum Number of Iterations have reached and ELBO has not yet converged\nThe difference in ELBO at stopping time is ", abs(diff_ELBO), "\nEBLO value at stop = ", elbo[iter], "\n"))
      Time.not.converged = Sys.time()
      Time.taken = round(Time.not.converged - Time.start, 3)
      Time.taken.unit = attr(Time.taken, "units")
      Status = "Not Converged"
    }
    # if( ((diff_ELBO<0) & (abs(diff_ELBO) > 1e-5))){
    #   cat(paste0("Warning! Iteration #", iter, " presents an ELBO decrement!\n"))
    # }
    if((abs(diff_ELBO) < epsilon) & (iter < maxIter)){
      # time at the convergence
      Time.converged = Sys.time()
      Time.taken = round(Time.converged - Time.start, 3)
      Time.taken.unit = attr(Time.taken, "units")
      Status = "Converged"
      cat(paste0("Success! Increase in ELBO < ", epsilon, ". Algorithm has converged in ",iter," iterations\n", "EBLO value at convergence = ", elbo[iter], "\nTime taken to converge: ", Time.taken," ", attr(Time.taken, "units")[1], "\n"))
      break
    }
    
  }
  
  est.s = apply(RHO_jk, 1, which.max)
  z.est = list()
  for(j in 1:J){
    z.est[[j]] = apply(XI_jil[[j]], 1, which.max)
  }
  
  z.estimated = NULL
  for(j in 1:J){
    z.estimated = c(z.estimated, z.est[[j]])
  }
  
  output = list(ELBO_val_conv = elbo[iter],
                ELBO = elbo,
                
                ml_y = ml_y,
                lambdal_y = lambdal_y,
                nul_y = nul_y,
                Sl_y = Sl_y,
                Wl_y = Wl_y,
                XI_jil = XI_jil,
                a_bar_Ulk = a_bar_Ulk,
                b_bar_Ulk = b_bar_Ulk,
                
                RHO_jk = RHO_jk,
                a_vk = a_vk,
                b_vk = b_vk,
                
                est.s = est.s,
                z.est = z.est,
                z.estimated = z.estimated,
                
                s1 = Update_alpha[1, ],
                s2 = Update_alpha[2, ],
                
                r1 = Update_beta[1,],
                r2 = Update_beta[2,],
                
                Time.taken = Time.taken,
                Time.taken.unit = Time.taken.unit,
                Status = Status
  )
  return(output)
}