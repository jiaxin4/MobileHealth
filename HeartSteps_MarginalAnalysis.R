# Note for reproducibility:
# To reproduce the data analysis results in the paper, please run the HeartSteps_Preprocessing.R first to get the data in 
# the form for analysis. To run the entire analysis (with cross-validation), we used a high performance cluster with multiple nodes and around
# approximately 70 core hours.



# functions needed --------------------------------------------------------
inCI <- function(x, upr, lwr) {
  all(x >= lwr & x <= upr)
}

rep_j <- function(mat, t_total.i = t_total.i){
  ncols <- ncol(mat)
  j_total.i <- nrow(mat)
  matrix(unlist(sapply(1:j_total.i, function(irow) rep(mat[irow,],t_total.i))),ncol = ncols, byrow= TRUE)
}

rep_t <- function(mat, j_total.i = j_total.i){
  do.call(rbind, replicate(j_total.i, mat, simplify=FALSE))
}


MatrixDiag <- function(M1,M2){
  rbind(cbind(M1, matrix(0,nrow=nrow(M1),ncol=ncol(M2))),
        cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
}

## "leave-one-out" K-fold CV 
CVdf_calculation <- function(dat, df1w, df2w, df1an, df2an, i_total, j_total_max, t_total_max, p1 = 0.3, p2 = 0.3){
  df1w.kloc <- attr(bs(1:j_total_max, df = df1w, intercept = FALSE),"knots")
  df2w.kloc <- attr(bs(1:t_total_max, df = df2w, intercept = FALSE),"knots")
  
  df1an.kloc <- attr(bs(1:j_total_max, df = df1an, intercept = FALSE),"knots")
  df2an.kloc <- attr(bs(1:t_total_max, df = df2an, intercept = FALSE),"knots")
  
  check <- c()
  cv <- 0
  user_list <- unique(dat$user)
  ##Use K-fold (K = 10) reduce computation time
  K <- 10
  #kfold_group <- sample(1:K, i_total, replace = TRUE, prob = rep(0.1, 10))
  list_hold <- c(rep(1:K, floor(i_total/K)), 1:(i_total-floor(i_total/K)*K))
  kfold_group <- sample(x = list_hold, size = i_total)
  kfold_group <- data.frame(user = user_list, kfold_group = kfold_group)
  dat <- merge(x = dat, y = kfold_group, by = "user")
  kgroup_num <- unique(kfold_group$kfold_group)
  for (groupi in kgroup_num){
    dat_remain <- dat %>% filter(kfold_group != groupi)
    dat_kgroup <- dat[dat$kfold_group == groupi, ]
    
    user_remain_list <- unique(dat_remain$user)
    gamma1 <- 0
    gamma2 <- 0
    for (i_ind in 1:length(user_remain_list)){
      user.i <- user_remain_list[i_ind]
      dat_i <- dat_remain[dat_remain$user == user.i,]
      
      j.vec.i <- unique(dat_i$decision.index.nogap)
      j_total.i <- length(j.vec.i)
      t.vec.i <- unique(dat_i$min.after.decision)
      t_total.i <- length(t.vec.i)
      bs.j.f1w.comb <- rep_j(bs(j.vec.i, knots = df1w.kloc, intercept = FALSE), t_total.i)
      bs.t.f2w.comb <- rep_t(bs(t.vec.i, knots = df2w.kloc, intercept = FALSE), j_total.i)
      
      bs.j.f1an.comb <- rep_j(bs(j.vec.i, knots = df1an.kloc, intercept = FALSE), t_total.i)
      bs.t.f2an.comb <- rep_t(bs(t.vec.i, knots = df2an.kloc, intercept = FALSE), j_total.i)
      
      U.i_w <- cbind(1, bs.j.f1w.comb, bs.t.f2w.comb)
      U.i_an <- cbind(1, bs.j.f1an.comb, bs.t.f2an.comb)
      U.i <- MatrixDiag(U.i_w, U.i_an)
      
      Atilde.i <- cbind(diag(dat_i$send.active - p1), diag(dat_i$send.sedentary - p2))
      gamma1.i <- t(U.i)%*%t(Atilde.i)%*%Atilde.i%*%U.i
      gamma1 <- gamma1 + gamma1.i
      
      y.i <- dat_i$y_minuseta
      gamma2.i <- t(U.i)%*%t(Atilde.i)%*%y.i
      gamma2 <- gamma2 + gamma2.i
    }
    gamma1_inv <- chol2inv(chol(gamma1))
    gamma_hat <-  gamma1_inv%*%gamma2
    fA1_spline <- c(rep("A1f(j)", times = df1w), rep("A1f(t)", times = df2w))
    fA2_spline <- c(rep("A2f(j)", times = df1an), rep("A2f(t)", times = df2an))
    split_gamma_hat <- split(gamma_hat, c("A1int", fA1_spline, "A2int", fA2_spline))
    
    #calculate the CV for all individuals in K-group
    kgroup_userlist <- unique(dat_kgroup$user)
    for (kgroup_i in kgroup_userlist){
      rem_data <- dat_kgroup[dat_kgroup$user == kgroup_i, ]
      j_total_remo <- length(unique(rem_data$decision.index.nogap))
      t_total_remo <- length(unique(rem_data$min.after.decision))
      
      bs.j.f1w <- rep_j(bs(1:j_total_remo, knots = df1w.kloc, intercept = FALSE), t_total_remo)
      bs.t.f2w <- rep_t(bs(1:t_total_remo, knots = df2w.kloc, intercept = FALSE), j_total_remo)
      
      bs.j.f1an <- rep_j(bs(1:j_total_remo, knots = df1an.kloc, intercept = FALSE), t_total_remo)
      bs.t.f2an <- rep_t(bs(1:t_total_remo, knots = df2an.kloc, intercept = FALSE), j_total_remo)
      
      bs_colA1 <- list(bs.j.f1w, bs.t.f2w)
      bs_colA2 <- list(bs.j.f1an, bs.t.f2an)
      names(bs_colA1) <- c("A1f(j)", "A1f(t)")
      names(bs_colA2) <- c("A2f(j)", "A2f(t)")
      bs_col <- c(bs_colA1, bs_colA2)
      f_hat_mat <- c()
      for (f_l in 1:length(bs_col)){
        f_name <- names(bs_col)[[f_l]]
        f_hat_mat <- cbind(f_hat_mat, bs_col[[f_name]] %*% split_gamma_hat[[f_name]])
        colnames(f_hat_mat)[f_l] <- f_name
      }
      f_hat_mat <- as.data.frame(f_hat_mat)
      
      cv_i <- sum((rem_data$y_minuseta - (rem_data$send.active - p1)*(split_gamma_hat$A1int + f_hat_mat$`A1f(j)` + f_hat_mat$`A1f(t)`)
                   - (rem_data$send.sedentary - p2)*(split_gamma_hat$A2int + f_hat_mat$`A2f(j)` + f_hat_mat$`A2f(t)`))^2)
      cv <- cv + cv_i
      check <- c(check, kgroup_i)
    }
  }
  output <- list(cv, list(df1w, df2w, df1an, df2an), check)
  output
}


# Select number of knots using K-fold cross validation--------------------------------------------------------
dta = readRDS("jbslot_public_60min.RDS")
library(tidyverse)
library(splines)
library(MASS)
require(lattice)
library(gridExtra)
library(mgcv)

dta$send.fac <- as.factor(ifelse(dta$send == TRUE, 1, 0))
dta$send.num <- as.numeric(ifelse(dta$send == TRUE, 1, 0))

dta$weekday <- weekdays(dta$decision.utime)
dta$is.weekday <- factor(dta$weekday %in% 
                           c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
dta$log.steps <- log(dta$steps + 0.5)

# Add availability variable
dta$avail.num <- as.numeric(ifelse(dta$avail == TRUE, 1, 0))

#construct some variables needed for estimating eta.
dta <- dta %>% 
  group_by(user) %>% 
  mutate(log.steps_lag1  = dplyr::lag(log.steps, 1, default = log(0.5)))

dat = dta

# Only take observations that J_i < 210
dat <- dat %>% filter(decision.index.nogap <= 210)
dat.dec.min <- dat %>% group_by(user) %>% summarise(max.dec = max(decision.index.nogap),
                                                    max.min = max(min.after.decision))
i_total <- length(unique(dat$user))
j_total_max <- 210 + 1
t_total_max <- min(dat.dec.min$max.min) + 1
p1 = p2 = 0.3

##estimate eta
fit_eta <- gam(log.steps ~ s(min.after.decision, bs = "cr") +                
                 s(decision.index.nogap, bs = "cr") +  
                 jbsteps30pre.log +  
                 location.other + 
                 is.weekday + 
                 log.steps_lag1, data = dat)
eta_hat <- predict(fit_eta, newdata = dat)
dat$eta_hat <- eta_hat
dat$y_minuseta <- dat$log.steps - dat$eta_hat

## calcualte CV ------------------------------------------------------------------------
df_design <- expand.grid(df1w = 3:10, df2w = 3:10, 
                         df1an = 3:10, df2an = 3:10)

cv_collector <- c()
for (iter in 1:dim(df_design)[1]){
  df1w <- df_design[iter, "df1w"]
  df2w <- df_design[iter, "df2w"]
  
  df1an <- df_design[iter, "df1an"]
  df2an <- df_design[iter, "df2an"]
  cv_collector[[iter]] <- CVdf_calculation(dat, df1w, df2w, df1an, df2an, i_total, j_total_max, t_total_max)
}

cv <- sapply(cv_collector, "[[", 1)
cv_dat <- as.data.frame(cbind(cv, t(sapply(cv_collector, "[[", 2))))
colnames(cv_dat) <- c("CV", "df1w", "df2w", 
                      "df1an", "df2an")
min_row <- which.min(cv_dat$CV)
selected_num_knots <- cv_dat[min_row, -1]

df1w_selected <- selected_num_knots$df1w
df2w_selected <- selected_num_knots$df2w
df1an_selected <- selected_num_knots$df1an
df2an_selected <- selected_num_knots$df2an


# Marginal Analysis with the selected number of knots--------------------------------------------------------

#number of knots
df1w = df1w_selected #6
df2w = df2w_selected #3
dfsum_w <- 1 + df1w + df2w
bs.j.f1w <- bs(1:j_total_max, df = df1w, intercept = FALSE)
bs.t.f2w <- bs(1:t_total_max, df = df2w, intercept = FALSE)
df1w.kloc <- attr(bs.j.f1w, "knots")
df2w.kloc <- attr(bs.t.f2w, "knots")

df1an = df1an_selected #3
df2an = df2an_selected #8
dfsum_an <- 1 + df1an + df2an
bs.j.f1an <- bs(1:j_total_max, df = df1an, intercept = FALSE)
bs.t.f2an <- bs(1:t_total_max, df = df2an, intercept = FALSE)
df1an.kloc <- attr(bs.j.f1an, "knots")
df2an.kloc <- attr(bs.t.f2an, "knots")

gamma1 <- 0
gamma2 <- 0
U_collect <- list()
user_list <- unique(dat$user)
for (i_ind in 1:length(user_list)){
  user.i <- user_list[i_ind]
  dat_i <- dat[dat$user == user.i,]
  
  j.vec.i <- unique(dat_i$decision.index.nogap)
  j_total.i <- length(j.vec.i)
  t.vec.i <- unique(dat_i$min.after.decision)
  t_total.i <- length(t.vec.i)
  bs.j.f1w.comb <- rep_j(bs(j.vec.i, knots = df1w.kloc, intercept = FALSE), t_total.i)
  bs.t.f2w.comb <- rep_t(bs(t.vec.i, knots = df2w.kloc, intercept = FALSE), j_total.i)
  
  bs.j.f1an.comb <- rep_j(bs(j.vec.i, knots = df1an.kloc, intercept = FALSE), t_total.i)
  bs.t.f2an.comb <- rep_t(bs(t.vec.i, knots = df2an.kloc, intercept = FALSE), j_total.i)
  
  U.i_w <- cbind(1, bs.j.f1w.comb, bs.t.f2w.comb)
  U.i_an <- cbind(1, bs.j.f1an.comb, bs.t.f2an.comb)
  U.i <- MatrixDiag(U.i_w, U.i_an)
  U_collect[[i_ind]] <- U.i
  
  #add availability variable
  I.i <- diag(dat_i$avail.num)
  
  Atilde.i <- cbind(diag(dat_i$send.active - p1), diag(dat_i$send.sedentary - p2))
  gamma1.i <- t(U.i)%*%t(Atilde.i)%*%I.i%*%Atilde.i%*%U.i
  gamma1 <- gamma1 + gamma1.i
  
  y.i <- dat_i$y_minuseta
  gamma2.i <- t(U.i)%*%t(Atilde.i)%*%I.i%*%y.i
  gamma2 <- gamma2 + gamma2.i
}
gamma1_inv <- chol2inv(chol(gamma1))
gamma_hat <-  gamma1_inv%*%gamma2
fA1_spline <- c(rep("A1f(j)", times = df1w), rep("A1f(t)", times = df2w))
fA2_spline <- c(rep("A2f(j)", times = df1an), rep("A2f(t)", times = df2an))
split_gamma_hat <- split(gamma_hat, c("A1int", fA1_spline, "A2int", fA2_spline))

bs_colA1 <- list(bs.j.f1w, bs.t.f2w)
bs_colA2 <- list(bs.j.f1an, bs.t.f2an)
names(bs_colA1) <- c("A1f(j)", "A1f(t)")
names(bs_colA2) <- c("A2f(j)", "A2f(t)")
bs_col <- c(bs_colA1, bs_colA2)
f_hat_mat <- c()
for (f_l in 1:length(bs_col)){
  f_name <- names(bs_col)[[f_l]]
  f_hat_mat[[f_l]] <- bs_col[[f_name]] %*% split_gamma_hat[[f_name]]
  names(f_hat_mat)[[f_l]] <- f_name
}
f_hat_mat <- append(f_hat_mat, split_gamma_hat$A1int)
f_hat_mat <- append(f_hat_mat, split_gamma_hat$A2int)
names(f_hat_mat) <- c(names(bs_col), "A1int", "A2int")

## Confidence intervals
asyvar_mid <- est_equ <- riri_corrected_all <- 0
for (i_ind in 1:length(user_list)){
  user.i <- user_list[i_ind]
  dat_i <- dat[dat$user == user.i,]
  Atilde.i <- cbind(diag(dat_i$send.active - p1), diag(dat_i$send.sedentary - p2))
  I.i <- diag(dat_i$avail.num)
  y.i <- dat_i$y_minuseta
  
  U.i <- U_collect[[i_ind]]
  
  est_equ.i <- t(U.i) %*% t(Atilde.i) %*% I.i %*% (y.i -  as.numeric(Atilde.i %*% U.i %*% gamma_hat))
  est_equ <- est_equ + est_equ.i
  asyvar_mid <- asyvar_mid + est_equ.i %*% t(est_equ.i)
  
  
  #small sample correction: calculate the middle of bias-corrected variance
  r.i <- I.i%*%(y.i - as.numeric(Atilde.i%*%U.i%*%gamma_hat))
  Hii <- I.i%*%Atilde.i%*%U.i%*%gamma1_inv%*%t(U.i)%*%t(Atilde.i)
  Ii <- diag(1, nrow = nrow(Hii), ncol = ncol(Hii))
  Ii_Hii_inv <- chol2inv(Ii - Hii)
  riri_corrected <- t(U.i)%*%t(Atilde.i)%*%Ii_Hii_inv%*%r.i%*%t(r.i)%*%t(Ii_Hii_inv)%*% Atilde.i%*%U.i
  riri_corrected_all <- riri_corrected_all + riri_corrected
}
bs_allA1 <- MatrixDiag(MatrixDiag(matrix(1),bs.j.f1w), bs.t.f2w)
bs_allA2 <- MatrixDiag(MatrixDiag(matrix(1),bs.j.f1an), bs.t.f2an)
bs_all <-  MatrixDiag(bs_allA1, bs_allA2)
asym_invmid <- gamma1_inv %*% asyvar_mid %*% t(gamma1_inv)
asyvar <- bs_all%*%asym_invmid%*%t(bs_all)
##bias-corrected variance
var_mid_corrected <- gamma1_inv%*%riri_corrected_all%*%t(gamma1_inv)
asycov_beta_corrected <- bs_all%*% var_mid_corrected %*%t(bs_all)


result <- list(list(gamma_hat, gamma1_inv), f_hat_mat, asym_invmid, bs_col, var_mid_corrected)
names(result) <- c("gamma_hat", "f_hat_mat", "asym_invmid", "bs_col", "var_mid_corrected")


# Make surface plots--------------------------------------------------------
dec.index.vec <- 1:j_total_max
mins.vec <- 1:t_total_max

## send.active surface plot
bs_col <- result$bs_col
bs.j.f1 <- bs_col$`A1f(j)`
bs.t.f2 <- bs_col$`A1f(t)`
int <- f_hat_mat$A1int
fj <- f_hat_mat$`A1f(j)`
ft <- f_hat_mat$`A1f(t)`

TrtEff <- TrtEff2 <- se <- c()
for (dec.index in dec.index.vec){
  TrtEff = c(TrtEff, int + fj[dec.index] + ft) #1
  
  bs_spec <- cbind(1, matrix(rep(bs.j.f1[dec.index,], each = t_total_max), ncol = df1w), 
                   bs.t.f2, 
                   matrix(0, nrow = t_total_max, ncol = dfsum_an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
print((max(abs(TrtEff - TrtEff2)) <= 10e-10))

mins_data <- data.frame(t.min = rep(1:t_total_max, times = length(dec.index.vec)))
mins_data$CEE <- TrtEff
mins_data$se <- se
mins_data$j.dec <- rep(dec.index.vec, each = t_total_max)

p1 <- plot_ly(mins_data, x = ~t.min, y = ~j.dec, z = ~CEE) %>% 
  add_markers(color = ~CEE) %>% 
  layout(title = paste0("Walking suggestion: "),
         scene = list(xaxis = list(title = "minute"), 
                      yaxis = list(title = "decision point"),
                      zaxis = list(range=c(-0.05, 0.15)))
  )
htmlwidgets::saveWidget(as_widget(p1), "surfacePlot_walking.html")


## Send.sedentary surface plot
int <- f_hat_mat$A2int
fj <- f_hat_mat$`A2f(j)`
ft <- f_hat_mat$`A2f(t)`
bs.j.f1 <- bs_col$`A2f(j)`
bs.t.f2 <- bs_col$`A2f(t)`

TrtEff <- TrtEff2 <- se <- c()
for (dec.index in dec.index.vec){
  TrtEff = c(TrtEff, int + fj[dec.index] + ft) #1
  
  bs_spec <- cbind(matrix(0, nrow = t_total_max, ncol = (dfsum_w)), 
                   1, matrix(rep(bs.j.f1[dec.index,], each = t_total_max), ncol = df1an), #2
                   bs.t.f2)
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
print((max(abs(TrtEff - TrtEff2)) <= 10e-10))

mins_data <- data.frame(t.min = rep(1:t_total_max, times = length(dec.index.vec)))
mins_data$CEE <- TrtEff
mins_data$se <- se
mins_data$j.dec <- rep(dec.index.vec, each = t_total_max)

p2 <- plot_ly(mins_data, x = ~t.min, y = ~j.dec, z = ~CEE) %>% 
  add_markers(color = ~CEE) %>% 
  layout(title = paste0("Anti-sedentary suggestion: "),
         scene = list(xaxis = list(title = "minute"), yaxis = list(title = "decision point")),
         yaxis = list(range=c(-0.05, 0.15)))
htmlwidgets::saveWidget(as_widget(p2), "surfacePlot_antisedentary.html")




































