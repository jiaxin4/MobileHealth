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


# "leave-one-out" K-fold CV 
CVdf_calculation <- function(dat, df1w, df2w, df3w, df4w, df1an, df2an, df3an, df4an, i_total, j_total_max, t_total_max, p1 = 0.3, p2 = 0.3){
  df1w.kloc <- attr(bs(1:j_total_max, df = df1w, intercept = FALSE),"knots")
  df2w.kloc <- attr(bs(1:t_total_max, df = df2w, intercept = FALSE),"knots")
  df3w.kloc <- attr(bs(1:j_total_max, df = df3w, intercept = FALSE),"knots")
  df4w.kloc <- attr(bs(1:t_total_max, df = df4w, intercept = FALSE),"knots")
  
  df1an.kloc <- attr(bs(1:j_total_max, df = df1an, intercept = FALSE),"knots")
  df2an.kloc <- attr(bs(1:t_total_max, df = df2an, intercept = FALSE),"knots")
  df3an.kloc <- attr(bs(1:j_total_max, df = df3an, intercept = FALSE),"knots")
  df4an.kloc <- attr(bs(1:t_total_max, df = df4an, intercept = FALSE),"knots")
  
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
      bs.j.f3w.comb <- rep_j(bs(j.vec.i, knots = df3w.kloc, intercept = FALSE), t_total.i)
      bs.t.f4w.comb <- rep_t(bs(t.vec.i, knots = df4w.kloc, intercept = FALSE), j_total.i)
      
      bs.j.f1an.comb <- rep_j(bs(j.vec.i, knots = df1an.kloc, intercept = FALSE), t_total.i)
      bs.t.f2an.comb <- rep_t(bs(t.vec.i, knots = df2an.kloc, intercept = FALSE), j_total.i)
      bs.j.f3an.comb <- rep_j(bs(j.vec.i, knots = df3an.kloc, intercept = FALSE), t_total.i)
      bs.t.f4an.comb <- rep_t(bs(t.vec.i, knots = df4an.kloc, intercept = FALSE), j_total.i)
      U.i_w <- cbind(1, bs.j.f1w.comb, bs.t.f2w.comb, 
                     bs.j.f3w.comb*dat_i$mod, bs.t.f4w.comb*dat_i$mod)
      U.i_an <- cbind(1, bs.j.f1an.comb, bs.t.f2an.comb, 
                      bs.j.f3an.comb*dat_i$mod, bs.t.f4an.comb*dat_i$mod)
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
    fA1_spline <- c(rep("A1f(j)", times = df1w), rep("A1f(t)", times = df2w), 
                    rep("A1f(j)*mod", times = df3w), rep("A1f(t)*mod", times = df4w))
    fA2_spline <- c(rep("A2f(j)", times = df1an), rep("A2f(t)", times = df2an), 
                    rep("A2f(j)*mod", times = df3an), rep("A2f(t)*mod", times = df4an))
    split_gamma_hat <- split(gamma_hat, c("A1int", fA1_spline, "A2int", fA2_spline))
    
    #calculate the CV for all individuals in K-group
    kgroup_userlist <- unique(dat_kgroup$user)
    for (kgroup_i in kgroup_userlist){
      rem_data <- dat_kgroup[dat_kgroup$user == kgroup_i, ]
      j_total_remo <- length(unique(rem_data$decision.index.nogap))
      t_total_remo <- length(unique(rem_data$min.after.decision))
      
      bs.j.f1w <- rep_j(bs(1:j_total_remo, knots = df1w.kloc, intercept = FALSE), t_total_remo)
      bs.t.f2w <- rep_t(bs(1:t_total_remo, knots = df2w.kloc, intercept = FALSE), j_total_remo)
      bs.j.f3w <- rep_j(bs(1:j_total_remo, knots = df3w.kloc, intercept = FALSE), t_total_remo)
      bs.t.f4w <- rep_t(bs(1:t_total_remo, knots = df4w.kloc, intercept = FALSE), j_total_remo)
      
      bs.j.f1an <- rep_j(bs(1:j_total_remo, knots = df1an.kloc, intercept = FALSE), t_total_remo)
      bs.t.f2an <- rep_t(bs(1:t_total_remo, knots = df2an.kloc, intercept = FALSE), j_total_remo)
      bs.j.f3an <- rep_j(bs(1:j_total_remo, knots = df3an.kloc, intercept = FALSE), t_total_remo)
      bs.t.f4an <- rep_t(bs(1:t_total_remo, knots = df4an.kloc, intercept = FALSE), j_total_remo)
      
      bs_colA1 <- list(bs.j.f1w, bs.t.f2w, bs.j.f3w, bs.t.f4w)
      bs_colA2 <- list(bs.j.f1an, bs.t.f2an, bs.j.f3an, bs.t.f4an)
      names(bs_colA1) <- c("A1f(j)", "A1f(t)", "A1f(j)*mod", "A1f(t)*mod")
      names(bs_colA2) <- c("A2f(j)", "A2f(t)", "A2f(j)*mod", "A2f(t)*mod")
      bs_col <- c(bs_colA1, bs_colA2)
      f_hat_mat <- c()
      for (f_l in 1:length(bs_col)){
        f_name <- names(bs_col)[[f_l]]
        f_hat_mat <- cbind(f_hat_mat, bs_col[[f_name]] %*% split_gamma_hat[[f_name]])
        colnames(f_hat_mat)[f_l] <- f_name
      }
      f_hat_mat <- as.data.frame(f_hat_mat)
      
      cv_i <- sum((rem_data$y_minuseta - (rem_data$send.active - p1)*(split_gamma_hat$A1int + f_hat_mat$`A1f(j)` + f_hat_mat$`A1f(t)` 
                                                                      + rem_data$mod*f_hat_mat$`A1f(j)*mod` + rem_data$mod*f_hat_mat$`A1f(t)*mod`)
                   - (rem_data$send.sedentary - p2)*(split_gamma_hat$A2int + f_hat_mat$`A2f(j)` + f_hat_mat$`A2f(t)` 
                                                     + rem_data$mod*f_hat_mat$`A2f(j)*mod` + rem_data$mod*f_hat_mat$`A2f(t)*mod`))^2)
      cv <- cv + cv_i
      check <- c(check, kgroup_i)
    }
  }
  output <- list(cv, list(df1w, df2w, df3w, df4w, df1an, df2an, df3an, df4an), check)
  output
}



# Select number of knots --------------------------------------------------------
dta = readRDS("jbslot_public_60min.RDS")
dta$send.fac <- as.factor(ifelse(dta$send == TRUE, 1, 0))
dta$send.num <- as.numeric(ifelse(dta$send == TRUE, 1, 0))

dta$weekday <- weekdays(dta$decision.utime)
dta$is.weekday <- factor(dta$weekday %in% 
                           c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
dta$log.steps <- log(dta$steps + 0.5)
dta$is.weekday.num <- ifelse(dta$is.weekday == TRUE, 1, 0)
dta$location.other.num <- ifelse(dta$location.other == TRUE, 1, 0)
dta$jbsteps30pre.zero.bin_56 <- ifelse(dta$jbsteps30pre.zero >= 56, 1, 0)

##The moderator for this analysis
dta$mod <- dta$location.other.num
plot_title_all <- c("other locations", "home/work")

# Add availability variable
dta$avail.num <- as.numeric(ifelse(dta$avail == TRUE, 1, 0))

#construct some variables needed for estimating eta.
dta <- dta %>% 
  group_by(user) %>% 
  mutate(log.steps_lag1  = dplyr::lag(log.steps, 1, default = log(0.5)))


## all individuals 
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
                 log.steps_lag1 +
                 jbsteps30pre.zero.bin_56, data = dat)
eta_hat <- predict(fit_eta, newdata = dat)
dat$eta_hat <- eta_hat
dat$y_minuseta <- dat$log.steps - dat$eta_hat

## calculate CV 
df_design <- expand.grid(df1w = 3:10, df2w = 3:10, df3w = 3:10, df4w = 3:10, 
                         df1an = 3:10, df2an = 3:10, df3an = 3:10, df4an = 3:10)

cv_collector <- c()
for (iter in 1:dim(df_design)[1]){
  df1w <- df_design[iter, "df1w"]
  df2w <- df_design[iter, "df2w"]
  df3w <- df_design[iter, "df3w"]
  df4w <- df_design[iter, "df4w"]
  
  df1an <- df_design[iter, "df1an"]
  df2an <- df_design[iter, "df2an"]
  df3an <- df_design[iter, "df3an"]
  df4an <- df_design[iter, "df4an"]
  cv_collector[[iter]] <- CVdf_calculation(dat, df1w, df2w, df3w, df4w, df1an, df2an, df3an, df4an, i_total, j_total_max, t_total_max)
}

## Find the knots number that has the smallest CV
cv <- sapply(cv_collector, "[[", 1)
cv_dat <- as.data.frame(cbind(cv, t(sapply(cv_collector, "[[", 2))))
colnames(cv_dat) <- c("CV", "df1w", "df2w", "df3w", "df4w", 
                      "df1an", "df2an", "df3an", "df4an")
min_row <- which.min(cv_dat$CV)
selected_num_knots <- cv_dat[min_row, -1]

df1w_selected <- selected_num_knots$df1w
df2w_selected <- selected_num_knots$df2w
df3w_selected <- selected_num_knots$df3w
df4w_selected <- selected_num_knots$df4w
df1an_selected <- selected_num_knots$df1an
df2an_selected <- selected_num_knots$df2an
df3an_selected <- selected_num_knots$df3an
df4an_selected <- selected_num_knots$df4an


# Data analysis with Moderators --------------------------------------------------------
df1w = df1w_selected #6
df2w = df2w_selected #4
df3w = df3w_selected #3
df4w = df4w_selected #10
dfsum_w <- 1 + df1w + df2w + df3w + df4w
bs.j.f1w <- bs(1:j_total_max, df = df1w, intercept = FALSE)
bs.t.f2w <- bs(1:t_total_max, df = df2w, intercept = FALSE)
bs.j.f3w <- bs(1:j_total_max, df = df3w, intercept = FALSE) 
bs.t.f4w <- bs(1:t_total_max, df = df4w, intercept = FALSE)
df1w.kloc <- attr(bs.j.f1w, "knots")
df2w.kloc <- attr(bs.t.f2w, "knots")
df3w.kloc <- attr(bs.j.f3w, "knots")
df4w.kloc <- attr(bs.t.f4w, "knots")

df1an = df1an_selected #3
df2an = df2an_selected #10
df3an = df3an_selected #6
df4an = df4an_selected #7
dfsum_an <- 1 + df1an + df2an + df3an + df4an
bs.j.f1an <- bs(1:j_total_max, df = df1an, intercept = FALSE)
bs.t.f2an <- bs(1:t_total_max, df = df2an, intercept = FALSE)
bs.j.f3an <- bs(1:j_total_max, df = df3an, intercept = FALSE) 
bs.t.f4an <- bs(1:t_total_max, df = df4an, intercept = FALSE)
df1an.kloc <- attr(bs.j.f1an, "knots")
df2an.kloc <- attr(bs.t.f2an, "knots")
df3an.kloc <- attr(bs.j.f3an, "knots")
df4an.kloc <- attr(bs.t.f4an, "knots")


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
  bs.j.f3w.comb <- rep_j(bs(j.vec.i, knots = df3w.kloc, intercept = FALSE), t_total.i)
  bs.t.f4w.comb <- rep_t(bs(t.vec.i, knots = df4w.kloc, intercept = FALSE), j_total.i)
  
  bs.j.f1an.comb <- rep_j(bs(j.vec.i, knots = df1an.kloc, intercept = FALSE), t_total.i)
  bs.t.f2an.comb <- rep_t(bs(t.vec.i, knots = df2an.kloc, intercept = FALSE), j_total.i)
  bs.j.f3an.comb <- rep_j(bs(j.vec.i, knots = df3an.kloc, intercept = FALSE), t_total.i)
  bs.t.f4an.comb <- rep_t(bs(t.vec.i, knots = df4an.kloc, intercept = FALSE), j_total.i)
 
  U.i_w <- cbind(1, dat_i$mod, bs.j.f1w.comb, bs.t.f2w.comb, 
               bs.j.f3w.comb*dat_i$mod, bs.t.f4w.comb*dat_i$mod)
  U.i_an <- cbind(1, dat_i$mod, bs.j.f1an.comb, bs.t.f2an.comb, 
                 bs.j.f3an.comb*dat_i$mod, bs.t.f4an.comb*dat_i$mod)
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
fA1_spline <- c(rep("A1f(j)", times = df1w), rep("A1f(t)", times = df2w), 
                rep("A1f(j)*mod", times = df3w), rep("A1f(t)*mod", times = df4w))
fA2_spline <- c(rep("A2f(j)", times = df1an), rep("A2f(t)", times = df2an), 
                rep("A2f(j)*mod", times = df3an), rep("A2f(t)*mod", times = df4an))
split_gamma_hat <- split(gamma_hat, c("A1int", "A1mod", fA1_spline, "A2int", "A2mod", fA2_spline))

bs_colA1 <- list(bs.j.f1w, bs.t.f2w, bs.j.f3w, bs.t.f4w)
bs_colA2 <- list(bs.j.f1an, bs.t.f2an, bs.j.f3an, bs.t.f4an)
names(bs_colA1) <- c("A1f(j)", "A1f(t)", "A1f(j)*mod", "A1f(t)*mod")
names(bs_colA2) <- c("A2f(j)", "A2f(t)", "A2f(j)*mod", "A2f(t)*mod")
bs_col <- c(bs_colA1, bs_colA2)
f_hat_mat <- c()
for (f_l in 1:length(bs_col)){
  f_name <- names(bs_col)[[f_l]]
  f_hat_mat[[f_l]] <- bs_col[[f_name]] %*% split_gamma_hat[[f_name]]
  names(f_hat_mat)[[f_l]] <- f_name
}
f_hat_mat <- c(f_hat_mat, split_gamma_hat$A1int, split_gamma_hat$A1mod,
               split_gamma_hat$A2int, split_gamma_hat$A2mod)
names(f_hat_mat) <- c(names(bs_col), "A1int", "A1mod", "A2int", "A2mod")

## Confidence intervals ----------------------------------------------------------------------
asyvar_mid <- est_equ <- riri_corrected_all <- 0
for (i_ind in 1:length(user_list)){
  user.i <- user_list[i_ind]
  dat_i <- dat[dat$user == user.i,]
  Atilde.i <- cbind(diag(dat_i$send.active - p1), diag(dat_i$send.sedentary - p2))
  I.i <- diag(dat_i$avail.num)
  y.i <- dat_i$y_minuseta
  
  U.i <- U_collect[[i_ind]]

  est_equ.i <- t(U.i) %*% t(Atilde.i) %*% I.i %*% (y.i -  as.numeric(Atilde.i %*% U.i %*% gamma_hat))
  #est_equ_collect[[i_ind]] <- est_equ.i
  est_equ <- est_equ + est_equ.i
  asyvar_mid <- asyvar_mid + est_equ.i %*% t(est_equ.i)
  
  
  #small sample correction: calculate the middle of bias-corrected variance
  r.i <- I.i%*%(y.i - as.numeric(Atilde.i%*%U.i%*%gamma_hat))
  Hii <- I.i%*%Atilde.i%*%U.i%*%gamma1_inv%*%t(U.i)%*%t(Atilde.i)
  Ii <- diag(1, nrow = nrow(Hii), ncol = ncol(Hii))
  Ii_Hii_inv <- chol2inv(Ii - Hii)
  riri_corrected <- t(U.i)%*%t(Atilde.i)%*%Ii_Hii_inv%*%r.i%*%t(r.i)%*%Ii_Hii_inv%*% Atilde.i%*%U.i
  riri_corrected_all <- riri_corrected_all + riri_corrected
}
bs_allA1 <- MatrixDiag(MatrixDiag(MatrixDiag(MatrixDiag(diag(c(1,1)),bs.j.f1w), bs.t.f2w), bs.j.f3w),bs.t.f4w)
bs_allA2 <- MatrixDiag(MatrixDiag(MatrixDiag(MatrixDiag(diag(c(1,1)),bs.j.f1an), bs.t.f2an), bs.j.f3an),bs.t.f4an)
bs_all <-  MatrixDiag(bs_allA1, bs_allA2)
asym_invmid <- gamma1_inv %*% asyvar_mid %*% t(gamma1_inv)
asyvar <- bs_all%*%asym_invmid%*%t(bs_all)
##bias-corrected variance
var_mid_corrected <- gamma1_inv%*%riri_corrected_all%*%t(gamma1_inv)
asycov_beta_corrected <- bs_all%*% var_mid_corrected %*%t(bs_all)


result <- list(list(gamma_hat, gamma1_inv), f_hat_mat, asym_invmid, bs_col, var_mid_corrected)
names(result) <- c("gamma_hat", "f_hat_mat", "asym_invmid", "bs_col", "var_mid_corrected")
#saveRDS(result, "result.RDS")


# Generate Plots --------------------------------------------------------
ylim.lower.min <- -0.2
ylim.upper.min <- 0.2
ylim.lower.dec <- -0.3
ylim.upper.dec <- 0.3
## Send.Walking ----------------------------------------------------------------------
asym_invmid <- var_mid_corrected #use corrected variance
bs_col <- result$bs_col
bs.j.f1 <- bs_col$`A1f(j)`
bs.t.f2 <- bs_col$`A1f(t)`
bs.j.f3 <- bs_col$`A1f(j)*mod`
bs.t.f4 <- bs_col$`A1f(t)*mod`
plot_mins <- c(1, 10, 15, 20, 40, 60)
plot_dec <- c(1, 20, 35, 40, 50, 70, 90, 180)
plot_days <- c(1, 4, 7, 8, 10, 14, 18, 36)

int <- f_hat_mat$A1int
slope_mod <- f_hat_mat$A1mod
fj <- f_hat_mat$`A1f(j)`
ft <- f_hat_mat$`A1f(t)`
fj_mod <- f_hat_mat$`A1f(j)*mod`
ft_mod <- f_hat_mat$`A1f(t)*mod`
### By mins, log.steps----------------------------------------------------------------------
library(RColorBrewer)
group_color <- brewer.pal(length(plot_dec), "Blues")
xlab_title <- "minute"
#### Mod = TRUE 
plot.title <- plot_title_all[1]
dec.index.vec <- plot_dec
TrtEff <- TrtEff2 <- se <- c()
for (dec.index in dec.index.vec){
  TrtEff = c(TrtEff, int + slope_mod + fj[dec.index] + ft + 
               fj_mod[dec.index] + ft_mod)
  
  bs_spec <- cbind(1, 1, matrix(rep(bs.j.f1[dec.index,], each = t_total_max), ncol = df1w), 
                   bs.t.f2, 
                   matrix(rep(bs.j.f3[dec.index,], each = t_total_max), ncol = df3w),
                   bs.t.f4,
                   matrix(0, nrow = t_total_max, ncol = dfsum_an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
mins_data <- data.frame(t.min = rep(1:t_total_max, times = length(dec.index.vec)))
mins_data$TrtEff <- TrtEff
mins_data$se <- se
mins_data$dec.index <- rep(dec.index.vec, each = t_total_max)

mins_data_dec1 <- mins_data
mins_data_dec1$moderator <- plot.title

### mod = FALSE 
plot.title <- plot_title_all[2]
TrtEff <- TrtEff2 <- se <- c()
for (dec.index in dec.index.vec){
  TrtEff = c(TrtEff, int + fj[dec.index] + ft) 
  
  bs_spec <- cbind(1, 0, matrix(rep(bs.j.f1[dec.index,], each = t_total_max), ncol = df1w), 
                   bs.t.f2, matrix(0, nrow = t_total_max, ncol = df3w),
                   matrix(0, nrow = t_total_max, ncol = df4w),
                   matrix(0, nrow = t_total_max, ncol = dfsum_an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
mins_data <- data.frame(t.min = rep(1:t_total_max, times = length(dec.index.vec)))
mins_data$TrtEff <- TrtEff
mins_data$se <- se
mins_data$dec.index <- rep(dec.index.vec, each = t_total_max)

mins_data_dec2 <- mins_data
mins_data_dec2$moderator <- plot.title

mins_data_dec_w <- list(mins_data_dec1, mins_data_dec2)

### By dec, log.steps ----------------------------------------------------------------------
group_color <-  brewer.pal(length(plot_mins), "YlGnBu")
xlab_title <- "decision point"
### mod = TRUE 
plot.title <- plot_title_all[1]
mins.vec <- plot_mins
TrtEff <- TrtEff2 <- se <- c()
for (mins in mins.vec){
  TrtEff = c(TrtEff, int + slope_mod + fj + ft[mins] + fj_mod + ft_mod[mins])
  
  bs_spec <- cbind(1, 1, bs.j.f1, 
                   matrix(rep(bs.t.f2[mins,], each = j_total_max), nrow = j_total_max),
                   bs.j.f3,
                   matrix(rep(bs.t.f4[mins,], each = j_total_max), nrow = j_total_max),
                   matrix(0, nrow = j_total_max, ncol = dfsum_an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
dec_data <- data.frame(j.dec = rep(1:j_total_max, times = length(mins.vec)))
dec_data$TrtEff <- TrtEff
dec_data$se <- se
dec_data$mins <- rep(mins.vec, each = j_total_max)

dec_data_min1 <- dec_data
dec_data_min1$moderator <- plot.title

### mod = FALSE 
plot.title <- plot_title_all[2]
mins.vec <- plot_mins
TrtEff <- TrtEff2 <- se <- c()
for (mins in mins.vec){
  TrtEff = c(TrtEff, int + fj + ft[mins])
  
  bs_spec <- cbind(1, 0, bs.j.f1, 
                   matrix(rep(bs.t.f2[mins,], each = j_total_max), nrow = j_total_max), 
                   matrix(0, ncol = df3w, nrow = j_total_max),
                   matrix(0, ncol = df4w, nrow = j_total_max),
                   matrix(0, nrow = j_total_max, ncol = dfsum_an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
dec_data <- data.frame(j.dec = rep(1:j_total_max, times = length(mins.vec)))
dec_data$TrtEff <- TrtEff
dec_data$se <- se
dec_data$mins <- rep(mins.vec, each = j_total_max)

dec_data_min2 <- dec_data
dec_data_min2$moderator <- plot.title

dec_data_min_w <- list(dec_data_min1, dec_data_min2)

## Send.sedentary ----------------------------------------------------------------------
int <- f_hat_mat$A2int
slope_mod <- f_hat_mat$A2mod
fj <- f_hat_mat$`A2f(j)`
ft <- f_hat_mat$`A2f(t)`
fj_mod <- f_hat_mat$`A2f(j)*mod`
ft_mod <- f_hat_mat$`A2f(t)*mod`
bs.j.f1 <- bs_col$`A2f(j)`
bs.t.f2 <- bs_col$`A2f(t)`
bs.j.f3 <- bs_col$`A2f(j)*mod`
bs.t.f4 <- bs_col$`A2f(t)*mod`

### By mins, log.steps----------------------------------------------------------------------
library(RColorBrewer)
group_color <- brewer.pal(length(plot_dec), "Blues")
xlab_title <- "minute"
#### mod = TRUE 
plot.title <- plot_title_all[1]
dec.index.vec <- plot_dec
TrtEff <- TrtEff2 <- se <- c()
for (dec.index in dec.index.vec){
  TrtEff = c(TrtEff, int + slope_mod + fj[dec.index] + ft + fj_mod[dec.index] + ft_mod) #1
  
  bs_spec <- cbind(matrix(0, nrow = t_total_max, ncol = (dfsum_w)), 
                   1, 1, matrix(rep(bs.j.f1[dec.index,], each = t_total_max), ncol = df1an), #2
                   bs.t.f2,
                   matrix(rep(bs.j.f3[dec.index,], each = t_total_max), ncol = df3an), 
                   bs.t.f4)
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
mins_data <- data.frame(t.min = rep(1:t_total_max, times = length(dec.index.vec)))
mins_data$TrtEff <- TrtEff
mins_data$se <- se
mins_data$dec.index <- rep(dec.index.vec, each = t_total_max)

mins_data_dec1 <- mins_data
mins_data_dec1$moderator <- plot.title

#### mod = FALSE 
plot.title <- plot_title_all[2]
TrtEff <- TrtEff2 <- se <- c()
for (dec.index in dec.index.vec){
  TrtEff = c(TrtEff, int + fj[dec.index] + ft) 
  
  bs_spec <- cbind(matrix(0, nrow = t_total_max, ncol = dfsum_w), 
                   1, 0, matrix(rep(bs.j.f1[dec.index,], each = t_total_max), ncol = df1an), 
                   bs.t.f2,
                   matrix(0, nrow = t_total_max, ncol = df3an),
                   matrix(0, nrow = t_total_max, ncol = df4an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
mins_data <- data.frame(t.min = rep(1:t_total_max, times = length(dec.index.vec)))
mins_data$TrtEff <- TrtEff
mins_data$se <- se
mins_data$dec.index <- rep(dec.index.vec, each = t_total_max)

mins_data_dec2 <- mins_data
mins_data_dec2$moderator <- plot.title

#### At each dec.index, four plots ----------------------------------------------------------------------
mins_data_dec_an <- list(mins_data_dec1, mins_data_dec2)
plot_list_w <- plot_list_an <- list()
require(lattice)
for (i_plot in 1:length(plot_dec)){
  dec.index.i <- plot_dec[i_plot]
  for (mod in 1:2){
    plot_dat <- mins_data_dec_w[[mod]] %>% filter(dec.index == dec.index.i)
    plot.title <- plot_dat$moderator[1]
    p_ci <- ggplot(plot_dat, aes(x = t.min, y = TrtEff)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_ribbon(aes(ymin = TrtEff - 1.96*se, ymax = TrtEff + 1.96*se), alpha = 0.3) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=11),
            axis.title=element_text(size=12),
            plot.subtitle=element_text(size=14)) +
      scale_color_manual(name = "Dec. index", values  = group_color[i_plot]) +
      scale_fill_manual(name = "Dec. index", values =  group_color[i_plot]) +
      coord_cartesian(ylim = c(ylim.lower.min, ylim.upper.min)) +
      labs(subtitle = plot.title) +
      xlab(xlab_title) +
      ylab("CEE")
    plot_list_w[[mod]] = p_ci
    
    plot_dat <- mins_data_dec_an[[mod]] %>% filter(dec.index == dec.index.i)
    plot.title <- plot_dat$moderator[1]
    p_ci <- ggplot(plot_dat, aes(x = t.min, y = TrtEff)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_ribbon(aes(ymin = TrtEff - 1.96*se, ymax = TrtEff + 1.96*se), alpha = 0.3) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=11),
            axis.title=element_text(size=12),
            plot.subtitle=element_text(size=14)) +
      scale_color_manual(name = "Dec. index", values  = group_color[i_plot]) +
      scale_fill_manual(name = "Dec. index", values =  group_color[i_plot]) +
      coord_cartesian(ylim = c(ylim.lower.min, ylim.upper.min)) +
      labs(subtitle = plot.title) +
      xlab(xlab_title) +
      ylab("CEE")
    plot_list_an[[mod]] = p_ci
  }
  p_comb_ci_w <- ggarrange(plot_list_w[[2]], plot_list_w[[1]], common.legend = TRUE)
  p_comb_ci_w <- annotate_figure(p_comb_ci_w,
                                 top=text_grob("Walking Suggestion", color = "black", face = "bold", size = 15))
  
  p_comb_ci_an <- ggarrange(plot_list_an[[2]], plot_list_an[[1]], common.legend = TRUE)
  p_comb_ci_an <- annotate_figure(p_comb_ci_an,
                                  top=text_grob("Anti-sedentary Suggestion", color = "black", face = "bold", size = 15))
  
  p_all <- ggarrange(p_comb_ci_w, p_comb_ci_an, ncol = 2, nrow = 1)
  pdf(file = paste0("CEE by mins at dec index ", dec.index.i," with CI.pdf"), width = 14, height = 3.5)
  print(p_all)
  dev.off()
}


### By dec, log.steps ----------------------------------------------------------------------
group_color <-  brewer.pal(length(plot_mins), "YlGnBu")
xlab_title <- "decision point"
#### mod = TRUE 
plot.title <- plot_title_all[1]
mins.vec <- plot_mins
TrtEff <- TrtEff2 <- se <- c()
for (mins in mins.vec){
  TrtEff = c(TrtEff, int + slope_mod + fj + ft[mins] + fj_mod + ft_mod[mins])
  
  bs_spec <- cbind(matrix(0, nrow = j_total_max, ncol = dfsum_w),
                   1, 1, bs.j.f1, 
                   matrix(rep(bs.t.f2[mins,], each = j_total_max), nrow = j_total_max), 
                   bs.j.f3,
                   matrix(rep(bs.t.f4[mins,], each = j_total_max), nrow = j_total_max))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
dec_data <- data.frame(j.dec = rep(1:j_total_max, times = length(mins.vec)))
dec_data$TrtEff <- TrtEff
dec_data$se <- se
dec_data$mins <- rep(mins.vec, each = j_total_max)

dec_data_min1 <- dec_data
dec_data_min1$moderator <- plot.title

#### mod = FALSE
plot.title <- plot_title_all[2]
mins.vec <- plot_mins
TrtEff <- TrtEff2 <- se <- c()
for (mins in mins.vec){
  TrtEff = c(TrtEff, int + fj + ft[mins])
  
  bs_spec <- cbind(matrix(0, nrow = j_total_max, ncol = dfsum_w),
                   1, 0, bs.j.f1, 
                   matrix(rep(bs.t.f2[mins,], each = j_total_max), nrow = j_total_max), 
                   matrix(0,ncol = df3an, nrow = j_total_max), 
                   matrix(0, nrow = j_total_max, ncol = df4an))
  TrtEff2 <- c(TrtEff2, bs_spec %*% gamma_hat) 
  
  se <- c(se, sqrt(diag(bs_spec%*%asym_invmid%*%t(bs_spec))))
}
dec_data <- data.frame(j.dec = rep(1:j_total_max, times = length(mins.vec)))
dec_data$TrtEff <- TrtEff
dec_data$se <- se
dec_data$mins <- rep(mins.vec, each = j_total_max)

dec_data_min2 <- dec_data
dec_data_min2$moderator <- plot.title

### At each min, four plots----------------------------------------------------------------------
dec_data_min_an <- list(dec_data_min1, dec_data_min2)
plot_list_w <- plot_list_an <- list()
require(lattice)
for (i_plot in 1:length(plot_mins)){
  mins.i <- plot_mins[i_plot]
  for (mod in 1:length(dec_data_min_an)){
    plot_dat <- dec_data_min_w[[mod]] %>% filter(mins == mins.i)
    plot.title <- plot_dat$moderator[1]
    p_ci <- ggplot(plot_dat, aes(x = j.dec, y = TrtEff)) +
      geom_line() +
      geom_ribbon(aes(ymin = TrtEff - 1.96*se, ymax = TrtEff + 1.96*se), alpha = 0.3) +
      theme_bw() +
      geom_hline(yintercept = 0, linetype = 3) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            plot.subtitle=element_text(size=14)) +
      scale_color_manual(name = "Mins", values  = group_color[i_plot]) +
      scale_fill_manual(name = "Mins", values =  group_color[i_plot]) +
      coord_cartesian(ylim = c(ylim.lower.dec, ylim.upper.dec)) +
      labs(subtitle = plot.title) +
      xlab(xlab_title) +
      ylab("CEE")
    plot_list_w[[mod]] = p_ci
    
    plot_dat <- dec_data_min_an[[mod]] %>% filter(mins == mins.i)
    plot.title <- plot_dat$moderator[1]
    p_ci <- ggplot(plot_dat, aes(x = j.dec, y = TrtEff)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_ribbon(aes(ymin = TrtEff - 1.96*se, ymax = TrtEff + 1.96*se), alpha = 0.3) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=11),
            axis.title=element_text(size=12),
            plot.subtitle=element_text(size=14)) +
      scale_color_manual(name = "Mins", values  = group_color[i_plot]) +
      scale_fill_manual(name = "Mins", values =  group_color[i_plot]) +
      coord_cartesian(ylim = c(ylim.lower.dec, ylim.upper.dec)) +
      labs(subtitle = plot.title) +
      xlab(xlab_title) +
      ylab("CEE")
    plot_list_an[[mod]] = p_ci
  }
  p_comb_ci_w <- ggarrange(plot_list_w[[2]], plot_list_w[[1]], common.legend = TRUE)
  p_comb_ci_w <- annotate_figure(p_comb_ci_w,
                                 top=text_grob("Walking Suggestion", color = "black", face = "bold", size = 15))
  
  p_comb_ci_an <- ggarrange(plot_list_an[[2]], plot_list_an[[1]], common.legend = TRUE)
  p_comb_ci_an <- annotate_figure(p_comb_ci_an,
                                  top=text_grob("Anti-sedentary Suggestion", color = "black", face = "bold", size = 15))
  
  p_all <- ggarrange(p_comb_ci_w, p_comb_ci_an, ncol = 2, nrow = 1)
  pdf(file = paste0("CEE by dec index at min ", mins.i," with CI.pdf"), width = 14, height = 3.5)
  print(p_all)
  dev.off()
}




















