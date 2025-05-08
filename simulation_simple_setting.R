# Note for reproducibility:
# To produce simulation study in the paper, we used a high performance cluster with multiple nodes and approximately 
# 200 core hours for 1000 simulations. 

#simulation
rm(list = ls())
source("simulation_function.R")
library(utils)
library(tidyverse)
library(splines)
library(MASS)
library(mgcv)
library(sn)

i_total<- 20 #in our actual simulation, we repeated the following process for other sample sizes, such as 50, 100.
nSeed <- 1 # We used 10 nodes at the same time and each node run 100 simulations, so the seed for each node is different. 
nsim <- 1000 
j_total = 30
t_total = 20

gen_formul <- "A*(alpha0_true + f1 + f2 + z1*f3 + b0_i_trt) + 0.35 + f1 + f2 + z1*g1 + g2*z3 + b0_i_nuis"

#function complexity
#complexity level: 1 < 2 < 3
f1_type <- "selfsimple"
f2_type <- "selfsimple"
f3_type <- 1


#settings
set.seed(nSeed)
pa = 0.4 #pmf of treatment
et_cor = 0 #correlation of e_i
alpha0_true <- 0.7

#generating functions
t_vec_eval <- seq(0.01,0.99,length = t_total)
j_vec_eval <- seq(0.01,0.99,length = j_total)
f1 <- gen_smooth_paper2(f1_type, j_vec_eval)
f2 <- gen_smooth_paper2(f2_type, t_vec_eval)
f3 <- gen_smooth_paper2(f3_type, t_vec_eval)
g1 <- gen_smooth_paper2(1, j_vec_eval)
g2 <- gen_smooth_paper2(3, t_vec_eval)


result_all <- list()
for (isim in 1:nsim){
  dta <- gen_data_discrete(i_total = i_total, j_total = j_total, t_total = t_total, gen_formul = gen_formul, 
                           f1_true = f1, f2_true = f2, f3_true = f3,
                           g1_true = g1, g2_true = g2,
                           et_cor = et_cor, pa = pa, alpha0_true = alpha0_true)
  
  # estimate eta  ----------------------------------------
  ##use GAM to fit.
  fit_eta <- gam(y ~ bs(j, df = 12) + bs(t, df = 12) + z1:bs(t, df = 12), data = dta)
  eta_hat <- predict(fit_eta, newdata = dta)
  dta$eta_hat <- eta_hat
  
  
  # Use cross-validation to select knots for each simulation
  dfselected <- eachsim_CVdf(dta = dta, i_total = i_total, j_total = j_total, t_total = t_total)
  df1_given <- as.numeric(dfselected$df1)
  df2_given <- as.numeric(dfselected$df2)
  df3_given <- as.numeric(dfselected$df3)
  
  # Estimate beta -----------------------------------------------------------
  ##B-spline but use parametric
  bs.j.f1 <- bs(j_vec_eval, df = df1_given, intercept = FALSE)
  bs.t.f2 <- bs(t_vec_eval, df = df2_given, intercept = FALSE)
  bs.t.f3 <- bs(t_vec_eval, df = df3_given, intercept = FALSE)
  bs.j.f1.comb <- rep_j(bs.j.f1, t_total)
  bs.t.f2.comb <- rep_t(bs.t.f2, j_total)
  bs.t.f3.comb <- rep_t(bs.t.f3, j_total)
  k_sum <- df1_given + df2_given + df3_given
  gamma1 <- matrix(0, ncol = k_sum + 1, nrow = k_sum + 1)
  gamma2 <- matrix(0, ncol = 1, nrow = k_sum + 1)
  U_collect <- list()
  for (i_ind in 1:i_total){
    dat_i <- dta[dta$id == i_ind,]
    U.i <- c()
    jt_total <- j_total*t_total
    
    U.i <- cbind(1, bs.j.f1.comb, bs.t.f2.comb, 
                 bs.t.f3.comb*dat_i$z1)
    
    U_collect[[i_ind]] <- U.i
    Atilde.i <- diag(dat_i$A_num - pa)
    gamma1.i <- t(U.i)%*%Atilde.i%*%Atilde.i%*%U.i
    gamma1 <- gamma1 + gamma1.i
    
    y.i <- dat_i$y - dat_i$eta_hat
    #here we set eta to be 0
    gamma2.i <- t(U.i)%*%Atilde.i%*%y.i
    gamma2 <- gamma2 + gamma2.i
  }
  
  gamma1_inv <- chol2inv(chol(gamma1))
  gamma_hat <-  gamma1_inv%*%gamma2
  #Verify the estimating equation = 0
  est_equ <- 0
  est_equ_collect <- c()
  for (i_ind in 1:i_total){
    dat_i <- dta[dta$id == i_ind,]
    Atilde.i <- diag(dat_i$A_num - pa)
    y.i <- dat_i$y - dat_i$eta_hat
    U.i <- U_collect[[i_ind]]
    est_equ.i <- t(U.i) %*% Atilde.i %*% (y.i -  as.numeric(Atilde.i %*% U.i %*% gamma_hat))
    est_equ_collect[[i_ind]] <- est_equ.i
    est_equ <- est_equ + est_equ.i
  }
  est_equ <- est_equ/i_total
  
  #f_hat
  alpha0_hat <- gamma_hat[1]
  f1_hat <- bs.j.f1%*%gamma_hat[2:(df1_given + 1)] 
  f2_hat <- bs.t.f2%*%gamma_hat[(df1_given+2):(df1_given + df2_given + 1)]
  f3_hat <- bs.t.f3%*%gamma_hat[(df1_given + df2_given + 2):(df1_given + df2_given + df3_given + 1)]
  #Calculate estimated variance
  beta_ce_hat <- alpha0_hat + rep(rep(f1_hat, each = t_total),times = i_total) + rep(rep(f2_hat, times = j_total), times = i_total) + dta$z1*rep(rep(f3_hat, times = j_total), times = i_total)
  epi <- (dta$A_num - pa)*(dta$y - dta$eta_hat) - pa*(1-pa)*beta_ce_hat
  ## Method 1: use all individuals' data to calculate V_i and use the same V_i for all i
  epi.matrix <- t(matrix(epi, nrow = t_total*j_total))
  cov_eta <- cov(epi.matrix)
  
  var_mid <- matrix(0, ncol = k_sum + 1, nrow = k_sum+1)
  gamma_tilde_1 <- matrix(0, ncol = k_sum+1, nrow = k_sum+1)
  gamma_tilde_2 <- matrix(0, ncol = 1, nrow = k_sum+1)
  gamma_hat_star_1 <- matrix(0, ncol = k_sum+1, nrow = k_sum+1)
  gamma_hat_star_2 <- matrix(0, ncol = 1, nrow = k_sum+1)
  var_mid_method2 <- 0
  riri_corrected_all <- 0
  for (i_ind in 1:i_total){
    dat_i <- dta[dta$id == i_ind,]
    U.i <- U_collect[[i_ind]]
    
    var_mid.i <- t(U.i)%*%cov_eta%*%U.i
    var_mid <- var_mid + var_mid.i
    
    beta_ce_i <- dat_i$f1 + dat_i$f2 + dat_i$z1*dat_i$f3
    gamma_tilde_1i <- pa*(1-pa)*t(U.i)%*%U.i
    gamma_tilde_1 <- gamma_tilde_1 + gamma_tilde_1i
    gamma_tilde_2i <- pa*(1-pa)*t(U.i)%*%beta_ce_i
    gamma_tilde_2 <- gamma_tilde_2 + gamma_tilde_2i
    
    Atilde.i <- diag(dat_i$A_num - pa)
    y.i <- dat_i$y - dat_i$eta_hat
    gamma_hat_star_2i <- t(U.i)%*%Atilde.i%*%y.i
    gamma_hat_star_2 <- gamma_hat_star_2 + gamma_hat_star_2i
    
    ## Method 2: sandwich estimator in TQ's notes
    est_equ.i <- est_equ_collect[[i_ind]]
    var_mid_method2 <- var_mid_method2 + est_equ.i %*% t(est_equ.i)
    
    ## Small sample correction
    r.i <- y.i - as.numeric(Atilde.i%*%U.i%*%gamma_hat)
    Hii <- Atilde.i%*%U.i%*%gamma1_inv%*%t(U.i)%*%Atilde.i
    Ii <- diag(1, nrow = nrow(Hii), ncol = ncol(Hii))
    Ii_Hii_inv <- chol2inv(Ii - Hii)
    riri_corrected <- t(U.i)%*%Atilde.i%*%Ii_Hii_inv%*%r.i%*%t(r.i)%*%Ii_Hii_inv%*% Atilde.i%*%U.i
    riri_corrected_all <- riri_corrected_all + riri_corrected
  }
  
  #f_tilde
  gamma_tilde_1_inv <- chol2inv(chol(gamma_tilde_1))
  gamma_tilde <- gamma_tilde_1_inv %*% gamma_tilde_2
  f1_tilde <- bs.j.f1%*%gamma_tilde[2:(df1_given + 1)] 
  f2_tilde <- bs.t.f2%*%gamma_tilde[(df1_given+2):(df1_given + df2_given + 1)]
  f3_tilde <- bs.t.f3%*%gamma_tilde[(df1_given + df2_given + 2):(df1_given + df2_given + df3_given + 1)]
  #f_hat_star
  gamma_hat_star <- gamma_tilde_1_inv %*% gamma_hat_star_2
  f1_star <- bs.j.f1%*%gamma_hat_star[2:(df1_given + 1)] 
  f2_star <- bs.t.f2%*%gamma_hat_star[(df1_given+2):(df1_given + df2_given + 1)]
  f3_star <- bs.t.f3%*%gamma_hat_star[(df1_given + df2_given + 2):(df1_given + df2_given + df3_given + 1)]
  
  #Variance
  ##1. Variance of beta_star
  bs_all <- MatrixDiag(MatrixDiag(MatrixDiag(matrix(1),bs.j.f1), bs.t.f2), bs.t.f3)
  cov_beta <- bs_all%*%gamma1_inv%*%var_mid%*%gamma1_inv%*%t(bs_all)
  cov_beta_star <- bs_all%*%gamma_tilde_1_inv%*%var_mid%*%gamma_tilde_1_inv%*%t(bs_all)
  ##2. Asymptotic variance
  asycov_beta <- bs_all%*%gamma_tilde_1_inv %*% var_mid_method2 %*% gamma_tilde_1_inv%*%t(bs_all)
  ###Asymptotic variance after small sample correction
  var_mid_corrected <- gamma1_inv%*%riri_corrected_all%*%gamma1_inv
  asycov_beta_corrected <- bs_all%*% var_mid_corrected %*%t(bs_all)
  
  
  #collect results of each f's
  f1_result <- evalcrit_summary(f = f1, fhat = f1_hat, fstar = f1_star, ftilde = f1_tilde,
                                asyvar_fhat = diag(asycov_beta)[2:(j_total+1)],
                                var_fhat = diag(asycov_beta_corrected)[2:(j_total+1)],
                                var_fstar = diag(cov_beta_star)[2:(j_total+1)])
  f2_result <- evalcrit_summary(f = f2, fhat = f2_hat, fstar = f2_star, ftilde= f2_tilde, 
                                asyvar_fhat =  diag(asycov_beta)[(j_total+2):(j_total + t_total+1)],
                                var_fhat = diag(asycov_beta_corrected)[(j_total+2):(j_total + t_total+1)], 
                                var_fstar = diag(cov_beta_star)[(j_total+2):(j_total + t_total+1)])
  f3_result <- evalcrit_summary(f = f3, fhat = f3_hat, fstar = f3_star, ftilde= f3_tilde, 
                                asyvar_fhat =  diag(asycov_beta)[(j_total + t_total + 2):(j_total + 2*t_total+1)],
                                var_fhat = diag(asycov_beta_corrected)[(j_total + t_total + 2):(j_total + 2*t_total+1)], 
                                var_fstar = diag(cov_beta_star)[(j_total + t_total + 2):(j_total + 2*t_total+1)])
  
  
  # Collect all results
  # Collect all results
  setting <- list(setting = data.frame(i_total = i_total, j_total = j_total, t_total = t_total, 
                                       df1 = df1_given, df2 = df2_given, df3 = df3_given), 
                  gen_formul = gen_formul, fresult_length = length(f3_result),
                  est_equ = est_equ)
  result <- c(f1_result, f2_result, f3_result, setting)
  
  
  #alpha0_hat result
  alpha0_var <- asycov_beta_corrected[1]
  alpha0_ci <- c(alpha0_hat - 1.96*sqrt(alpha0_var), alpha0_hat + 1.96*sqrt(alpha0_var))
  alpha0_in <- ifelse(alpha0_true >= alpha0_ci[1] & alpha0_true <= alpha0_ci[2], 1, 0)
  alpha0_result <- list(alpha0_hat, alpha0_var, alpha0_ci, alpha0_in)
  names(alpha0_result) <- c("alpha0_hat", "alpha0_var", "alpha0_ci", "alpha0_in")
  
  result <- append(result, alpha0_result)
  result_all[[isim]] <- result
}

saveRDS(result_all, paste0("simresult_", i_total, "_", nSeed, ".RDS"))









