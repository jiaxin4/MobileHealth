#functions
inCI <- function(x, upr, lwr) {
  all(x >= lwr & x <= upr)
}


# functions needed --------------------------------------------------------

#This function is same as gen_smooth_paper except that it takes the time_vec as input instead of time_total
gen_smooth_paper2 <- function(case_num, time_vec){
  if (case_num == 1){
    g <- (1/3)*dbeta(time_vec, 10,5) + (1/3)*dbeta(time_vec, 7,7) + (1/3)*dbeta(time_vec, 5,10) 
  } else if (case_num == 2){
    g <- (6/10)*dbeta(time_vec, 30,17) + (4/10)*dbeta(time_vec, 3,11) 
  } else if (case_num == 3){
    g <- (1/3)*dbeta(time_vec, 20,5) + (1/3)*dbeta(time_vec, 12, 12) + (1/3)*dbeta(time_vec, 7,30) 
  } else if (case_num == "sin") {
    g <- (sin(2*pi*(time_vec - 0.5)))^2
  } else if (case_num == "selfsimple"){
    g <- (1/3)*dbeta(time_vec, 8,5) + (1/3)*dbeta(time_vec, 7,7) + (1/3)*dbeta(time_vec, 5,8)
  } 
  #g <- g - mean(g) #identifiability issue
  g <- g - g[1]
}



gen_data_discrete <- function(i_total, t_total, j_total, gen_formul, 
                              f1_true, f2_true, f3_true, g1_true, g2_true,
                              et_cor = 0.4, pa = 0.4, alpha0_true){
  dta <- c()
  #yj_cov <- matrix(ej_cor, ncol = j_total, nrow = j_total) + diag(-ej_cor + 1, ncol = j_total, nrow= j_total)
  yt_cov <- matrix(et_cor, ncol = t_total, nrow = t_total) + diag(-et_cor + 1, ncol = t_total, nrow= t_total)
  for (i in 1:i_total){
    t.vec = seq(0.01, 0.99, length = t_total)
    j.vec = seq(0.01, 0.99, length = j_total)
    f1.i <- f1_true
    f2.i <- f2_true
    f3.i <- f3_true
    
    g1.i <- g1_true
    g2.i <- g2_true
    #1. fixed j, assign excangeable for all time point; 
    #2. treat t and j as the same
    A_ij <- rbinom(j_total, size = 1, p = pa)
    e_it <- MASS::mvrnorm(j_total, mu = rep(0, t_total), Sigma = yt_cov)
    
    z1_i <- rbinom(j_total, size = 1, p = 0.2)
    z2_i <- 1
    z3_i <- rnorm(j_total)
    
    dta.i <- data.frame(id = i, j = rep(j.vec, each = t_total), t = rep(t.vec, times = j_total), 
                        A_num = rep(A_ij, each = t_total), A = as.factor(rep(A_ij, each = t_total)),
                        f1 = rep(f1.i, each = t_total), f2 = rep(f2.i, times = j_total), f3 = rep(f3.i, times = j_total),
                        g1 = rep(g1.i, each = t_total), g2 = rep(g2.i, times = j_total),
                        z1 = rep(z1_i, times = t_total), z2 = rep(z2_i, times = t_total), z3 = rep(z3_i, times = t_total),
                        e = as.vector(t(e_it)))
    dta <- rbind(dta, dta.i)
  }
  
  ########################  
  A <- dta$A_num
  f1 <- dta$f1 
  f2 <- dta$f2
  f3 <- dta$f3
  g1 <- dta$g1
  g2 <- dta$g2
  z1 <- dta$z1
  z2 <- dta$z2
  z3 <- dta$z3
  e <- dta$e
  expect_y <- eval(parse(text=gen_formul))
  y <- expect_y + e
  #######################
  dta <- cbind(dta,expect_y, y)
  dta
}


rep_j <- function(mat, t_total = t_total){
  ncols <- ncol(mat)
  matrix(unlist(sapply(1:j_total, function(irow) rep(mat[irow,],t_total))),ncol = ncols, byrow= TRUE)
}

rep_t <- function(mat, j_total = j_total){
  do.call(rbind, replicate(j_total, mat, simplify=FALSE))
}

MatrixDiag <- function(M1,M2){
  rbind(cbind(M1, matrix(0,nrow=nrow(M1),ncol=ncol(M2))),
        cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
}


cp_general <- function(fin, f_est, var_f_ci){
  eval_points <- length(fin)
  ci_f_lower <- f_est - 1.96*sqrt(var_f_ci)
  ci_f_upper <- f_est + 1.96*sqrt(var_f_ci)
  ci_f <- list(ci_f_lower, ci_f_upper)
  f_in_star <-  ifelse(sapply(1:eval_points, function(p) any(ci_f_lower[p] <= fin[p] & ci_f_upper[p] >= fin[p])),1, 0)
  f_in_star
}

f_cp_summary <- function(f, fhat, fstar, ftilde, var_fstar,var_fhat, asyvar_fhat){
  fin_fhat_asyvarfhat <- cp_general(fin = f, f_est = fhat, var_f_ci = asyvar_fhat)
  fin_fhat_varfhat <- cp_general(fin = f, f_est = fhat, var_f_ci = var_fhat)
  fin_fhat_varfstar <- cp_general(fin = f, f_est = fhat, var_f_ci = var_fstar)
  fin_fstar_varfstar <- cp_general(fin = f, f_est = fstar, var_f_ci = var_fstar)
  ftildein_fstar_varfstar <- cp_general(fin = ftilde, f_est = fstar, var_f_ci = var_fstar)
  list(fin_fhat_asyvarfhat, fin_fhat_varfstar, fin_fhat_varfhat, fin_fstar_varfstar, ftildein_fstar_varfstar)
}

evalcrit_summary <- function(f, fhat, fstar, ftilde, asyvar_fhat, var_fhat, var_fstar){
  #MISE
  ## {f - f_hat}
  mise_l1_f_fhat <- sum(abs(f-fhat))
  mise_l2_f_fhat <- sqrt(sum((f - fhat)^2))
  ## {fhat - fstar}
  mise_l1_fhat_fstar <- sum(abs(fhat - fstar))
  mise_l2_fhat_fstar <- sqrt(sum((fhat - fstar)^2))
  ## {f_tilde - f_hat}
  mise_l1_ftilde_fhat <- sum(abs(fhat-ftilde))
  mise_l2_ftilde_fhat <- sqrt(sum((fhat - ftilde)^2))
  ## {f - f_tilde}
  mise_l1_ftilde_f <- sum(abs(f - ftilde))
  mise_l2_ftilde_f <- sqrt(sum((f - ftilde)^2))
  mise <- c(mise_l1_f_fhat, mise_l2_f_fhat,
            mise_l1_fhat_fstar, mise_l2_fhat_fstar,
            mise_l1_ftilde_fhat, mise_l2_ftilde_fhat,
            mise_l1_ftilde_f, mise_l2_ftilde_f)
  
  #Avergae Bias
  avgbias_f_fhat <- sum(fhat - f)
  avgbias_f_ftilde <- sum(ftilde - f)
  avgbias_fhat_ftilde <- sum(fhat - ftilde)
  avgbias <- c(avgbias_f_fhat, avgbias_f_ftilde, avgbias_fhat_ftilde)
  
  #pointwise convergence
  ptconv_f_fhat <- abs(f - fhat)
  ptconv_fstar_ftilde <- abs(fstar - ftilde)
  ptconv_fhat_ftilde <- abs(fhat - ftilde)
  ptconv_fhat_fstar <- abs(fhat - fstar)
  
  #bias
  bias_fhat_f <- fhat - f
  bias_fhat_fstar <- fhat - fstar
  bias_fhat_ftilde <- fhat - ftilde
  
  #CP
  f_cp <-  f_cp_summary(f = f, fhat = fhat, fstar = fstar, ftilde= ftilde, 
                        asyvar_fhat = asyvar_fhat,
                        var_fhat = var_fhat, 
                        var_fstar = var_fstar)
  
  #collect all result
  hold_list <- list(f, fhat,fstar, ftilde, 
                    mise, avgbias, 
                    ptconv_f_fhat, ptconv_fstar_ftilde, ptconv_fhat_ftilde, ptconv_fhat_fstar,
                    bias_fhat_f, bias_fhat_fstar, bias_fhat_ftilde,
                    asyvar_fhat, var_fhat, var_fstar)
  result_list <- append(f_cp, hold_list)
  result_list
}


output_organize <- function(result_all, start_index, vec_eval){
  fin_fhat_asyvarfhat <- sapply(result_all, "[[", start_index)
  fin_fhat_varfstar <- sapply(result_all, "[[", start_index+1)
  fin_fhat_varfhat <- sapply(result_all, "[[", start_index+2)
  fin_fstar_varfstar <- sapply(result_all, "[[", start_index+3)
  ftildein_fstar_varfstar <- sapply(result_all, "[[", start_index+4)
  f_all <- sapply(result_all, "[[", start_index+5)
  fhat_all <- sapply(result_all, "[[", start_index+6)
  fstar_all <- sapply(result_all, "[[", start_index+7)
  ftilde_all <- sapply(result_all, "[[", start_index+8)
  
  #output1: f CP table
  fCP_fhat_asyvarfhat <- apply(fin_fhat_asyvarfhat, 1, mean)
  fCP_fhat_varfstar <-  apply(fin_fhat_varfstar, 1, mean)
  fCP_fhat_varfhat <- apply(fin_fhat_varfhat, 1, mean)
  fCP_fstar_varfstar <- apply(fin_fstar_varfstar, 1, mean)
  ftildeCP_fstar_varfstar <- apply(ftildein_fstar_varfstar, 1, mean)
  ###calculate sample variance
  samvar_fhat <- apply(fhat_all, 1, var)
  f <- f_all[,1]
  ci_f_lower_all <- fhat_all - 1.96*sqrt(samvar_fhat)
  ci_f_upper_all <- fhat_all + 1.96*sqrt(samvar_fhat)
  samsd_f_in_all <- c()
  nsim <- dim(fhat_all)[2]
  for (icol in 1:nsim){
    samsd_f_in <-  ifelse(sapply(1:length(f), function(p) any(ci_f_lower_all[p, icol] <= f[p] & ci_f_upper_all[p, icol] >= f[p])),1, 0)
    samsd_f_in_all <- cbind(samsd_f_in_all, samsd_f_in)
  }
  fCP_fhat_samsd <- apply(samsd_f_in_all, 1, mean)
  ###all CP
  cp_dat <- data.frame(cbind(fCP_fhat_asyvarfhat, fCP_fhat_varfstar, fCP_fhat_varfhat, 
                             fCP_fstar_varfstar, ftildeCP_fstar_varfstar,
                             fCP_fhat_samsd, vec_eval))
  
  #output2: f and f_hat plots
  plot_dat_f <- c()
  for (k in 1:9){
    hold_dat_f <- data.frame(f, f_hat = fhat_all[,k], 
                             time = vec_eval, iter = rep(k, times = length(f)))
    plot_dat_f <- rbind(hold_dat_f, plot_dat_f)
  }
  
  #output3: setting
  setting <- result_all[[1]]$setting
  gen_formul <- result_all[[1]]$gen_formul
  
  #output4: mise, avgbias, avgCP
  mise_all <- apply(sapply(result_all, "[[", start_index + 9), 1, mean)
  avgbias_all <- apply(sapply(result_all, "[[", start_index + 10), 1,mean)
  avgCP_all <- mean(fCP_fhat_asyvarfhat)
  
  #output5: pointwise convergence
  ptconv_f_fhat <- apply(sapply(result_all, "[[", start_index + 11), 1, mean)
  ptconv_fstar_ftilde <- apply(sapply(result_all, "[[", start_index + 12), 1, mean)
  ptconv_fhat_ftilde <- apply(sapply(result_all, "[[", start_index + 13), 1, mean)
  ptconv_fhat_fstar <- apply(sapply(result_all, "[[", start_index + 14), 1, mean)
  ptconv_mse_f_fhat <- apply((f_all - fhat_all)^2, 1, mean)
  ptconv_dat <- cbind(f_fhat = ptconv_f_fhat, 
                      fstar_ftilde = ptconv_fstar_ftilde, 
                      fhat_ftilde = ptconv_fhat_ftilde, 
                      fhat_fstar = ptconv_fhat_fstar,
                      mse_f_fhat = ptconv_mse_f_fhat)
  
  #output6: bias
  bias_fhat_f <- apply(sapply(result_all, "[[", start_index+ 15), 1, mean)
  bias_fhat_fstar <- apply(sapply(result_all, "[[",start_index +16), 1, mean)
  bias_fhat_ftilde <- apply(sapply(result_all, "[[", start_index+17), 1, mean)
  bias_dat <- cbind(fhat_f = bias_fhat_f, 
                    fhat_fstar = bias_fhat_fstar, 
                    fhat_ftilde = bias_fhat_ftilde)
  
  #output7: variance
  asyvar_fhat <- apply(sapply(result_all, "[[", start_index+18), 1, mean)
  var_fhat <- apply(sapply(result_all, "[[", start_index+19), 1, mean)
  var_fstar <- apply(sapply(result_all, "[[", start_index+20), 1, mean)
  samvar_fstar <- apply(fstar_all, 1, var)
  samvar_fstar_ftilde <- apply(fstar_all - ftilde_all, 1, var)
  samvar_fhat_ftilde <- apply(fhat_all - ftilde_all, 1, var)
  var_dat <- cbind(fhat_asyvar = asyvar_fhat, 
                   fhat = var_fhat,
                   fstar = var_fstar,
                   fhat_samvar = samvar_fhat, 
                   fstar_samvar = samvar_fstar,
                   fstar_ftilde = samvar_fstar_ftilde,
                   fhat_ftilde = samvar_fhat_ftilde)
  
  result_f <- list(cp_dat = cp_dat, plot_dat_f = plot_dat_f, setting = setting, gen_formul = gen_formul, 
                   mise = mise_all, avgbias = avgbias_all, avgCP = avgCP_all,
                   ptconv = ptconv_dat, bias = bias_dat, var = var_dat)
  result_f
}


# "leave-one-out" K-fold CV --------------------------------------------------------
CVdf_calculation <- function(dat, df1, df2, df3, i_total, j_total, t_total, pa){
  t_vec_eval <- seq(0.01,0.99,length = t_total)
  j_vec_eval <- seq(0.01,0.99,length = j_total)
  df1.kloc <- attr(bs(j_vec_eval, df = df1, intercept = FALSE),"knots")
  df2.kloc <- attr(bs(t_vec_eval, df = df2, intercept = FALSE),"knots")
  df3.kloc <- attr(bs(j_vec_eval, df = df3, intercept = FALSE),"knots")
  
  check <- c()
  cv <- 0
  user_list <- unique(dat$id)
  ##Use K-fold (K = 10) reduce computation time
  K <- 10
  #kfold_group <- sample(1:K, i_total, replace = TRUE, prob = rep(0.1, 10))
  if (ceiling(i_total/K) == floor(i_total/K)){
    list_hold <- rep(1:K, floor(i_total/K))
  } else {
    list_hold <- c(rep(1:K, floor(i_total/K)), 1:(i_total-floor(i_total/K)*K))
  }
  kfold_group <- sample(x = list_hold, size = i_total)
  kfold_group <- data.frame(id = user_list, kfold_group = kfold_group)
  dat <- merge(x = dat, y = kfold_group, by = "id")
  kgroup_num <- unique(kfold_group$kfold_group)
  
  for (groupi in kgroup_num){
    dat_remain <- dat %>% filter(kfold_group != groupi)
    dat_kgroup <- dat[dat$kfold_group == groupi, ]
    
    user_remain_list <- unique(dat_remain$id)
    gamma1 <- 0
    gamma2 <- 0
    for (i_ind in 1:length(user_remain_list)){
      user.i <- user_remain_list[i_ind]
      dat_i <- dat_remain[dat_remain$id == user.i,]
      
      j.vec.i <- unique(dat_i$j)
      j_total.i <- length(j.vec.i)
      t.vec.i <- unique(dat_i$t)
      t_total.i <- length(t.vec.i)
      bs.j.f1.comb <- rep_j(bs(j.vec.i, knots = df1.kloc, intercept = FALSE), t_total.i)
      bs.t.f2.comb <- rep_t(bs(t.vec.i, knots = df2.kloc, intercept = FALSE), j_total.i)
      bs.t.f3.comb <- rep_t(bs(t.vec.i, knots = df3.kloc, intercept = FALSE), j_total.i)
      U.i <- cbind(1, bs.j.f1.comb, bs.t.f2.comb, 
                   bs.t.f3.comb*dat_i$z1)
      
      Atilde.i <- diag(dat_i$A_num - pa)
      gamma1.i <- t(U.i)%*%Atilde.i%*%Atilde.i%*%U.i
      gamma1 <- gamma1 + gamma1.i
      
      y.i <- dat_i$y - dat_i$eta_hat
      gamma2.i <- t(U.i)%*%Atilde.i%*%y.i
      gamma2 <- gamma2 + gamma2.i
    }
    gamma1_inv <- chol2inv(chol(gamma1))
    gamma_hat <-  gamma1_inv%*%gamma2
    fA_spline <- c(rep("Af(j)", times = df1), rep("Af(t)", times = df2),
                   rep("Af(t)*z1", times = df3))
    split_gamma_hat <- split(gamma_hat, c("Aint", fA_spline))
    
    #calculate the CV for all individuals in kth group
    kgroup_userlist <- unique(dat_kgroup$id)
    for (kgroup_i in kgroup_userlist){
      rem_data <- dat_kgroup[dat_kgroup$id == kgroup_i, ]

      bs_col <-  list(bs.j.f1.comb, bs.t.f2.comb, bs.t.f3.comb)
      names(bs_col) <- c("Af(j)", "Af(t)", "Af(t)*z1")
      f_hat_mat <- c()
      for (f_l in 1:length(bs_col)){
        f_name <- names(bs_col)[[f_l]]
        f_hat_mat <- cbind(f_hat_mat, bs_col[[f_name]] %*% split_gamma_hat[[f_name]])
        colnames(f_hat_mat)[f_l] <- f_name
      }
      f_hat_mat <- as.data.frame(f_hat_mat)
      
      cv_i <- sum((rem_data$y - rem_data$eta_hat - (rem_data$A_num - pa)*(split_gamma_hat$Aint + f_hat_mat$`Af(j)` + f_hat_mat$`Af(t)` + rem_data$z1*f_hat_mat$`Af(t)*z1`))^2)
      cv <- cv + cv_i
      check <- c(check, kgroup_i)
    }
  }
  output <- list(cv, list(df1, df2, df3), check)
  output
}

# for each simulation, calculate the cross-validation knots --------------------------------------------------------
eachsim_CVdf <- function(dta, i_total, j_total, t_total, pa = 0.4){
  # select df1 --------------------------------------------------------------
  df_given <- 7:18
  df_design <- expand.grid(df1 = df_given, df2 = 10, df3 = 10)
  cv_collector <- c()
  for (iter in 1:dim(df_design)[1]){
    df1 <- df_design[iter, "df1"]
    df2 <- df_design[iter, "df2"]
    df3 <- df_design[iter, "df3"]
    cv_collector[[iter]] <- CVdf_calculation(dat = dta, df1, df2, df3, i_total, j_total, t_total, pa = pa)
  }
  
  ## check smallest CV
  cv <- sapply(cv_collector, "[[", 1)
  cv_dat <- as.data.frame(cbind(cv, t(sapply(cv_collector, "[[", 2))))
  colnames(cv_dat) <- c("CV", "df1", "df2", "df3")
  dfselected <- cv_dat[which(cv_dat$CV == min(as.numeric(cv_dat$CV))), ]
  
  # select df2 --------------------------------------------------------------
  df1 <- as.numeric(dfselected$df1)
  df_design <- expand.grid(df1 = df1, df2 = df_given, df3 = 10)
  cv_collector <- c()
  for (iter in 1:dim(df_design)[1]){
    df1 <- df_design[iter, "df1"]
    df2 <- df_design[iter, "df2"]
    df3 <- df_design[iter, "df3"]
    cv_collector[[iter]] <- CVdf_calculation(dat = dta, df1, df2, df3, i_total, j_total, t_total, pa = pa)
  }
  
  ## check smallest CV
  cv <- sapply(cv_collector, "[[", 1)
  cv_dat <- as.data.frame(cbind(cv, t(sapply(cv_collector, "[[", 2))))
  colnames(cv_dat) <- c("CV", "df1", "df2", "df3")
  dfselected <- cv_dat[which(cv_dat$CV == min(as.numeric(cv_dat$CV))), ]
  
  # select df3 --------------------------------------------------------------
  df2 <- as.numeric(dfselected$df2)
  df_given <- 7:18
  df_design <- expand.grid(df1 = df1, df2 = df2, df3 = df_given)
  cv_collector <- c()
  for (iter in 1:dim(df_design)[1]){
    df1 <- df_design[iter, "df1"]
    df2 <- df_design[iter, "df2"]
    df3 <- df_design[iter, "df3"]
    cv_collector[[iter]] <- CVdf_calculation(dat = dta, df1, df2, df3, i_total, j_total, t_total, pa = pa)
  }
  
  ## check smallest CV
  cv <- sapply(cv_collector, "[[", 1)
  cv_dat <- as.data.frame(cbind(cv, t(sapply(cv_collector, "[[", 2))))
  colnames(cv_dat) <- c("CV", "df1", "df2", "df3")
  dfselected <- cv_dat[which(cv_dat$CV == min(as.numeric(cv_dat$CV))), ]
  
  dfselected
}










