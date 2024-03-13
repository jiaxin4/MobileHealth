rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(knitr)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(latex2exp)

simu_design <- data.frame(i_total = c(20, 50, 100),
                          j_total = 30,
                          t_total = 20)
seedlist <- list(1:10, 1:10, 1:10)
source("sim_function_CV.R")


# Organize results -------------------------------------------

result_all <- c()
for (i_design in 1:dim(simu_design)[1]){
  i_total <- simu_design[i_design,]$i_total
  seedlist.i <- seedlist[[i_design]]
  result_all <- c()
  for (iSeed in seedlist.i){
    nSeed <- iSeed
    result <- readRDS(paste0("simresult_", i_total, "_", nSeed, ".RDS"))
    result_all <- append(result_all, result)
  }
  print(length(result_all))
  fresult_length <- result_all[[1]]$fresult_length
  t_vec_eval <- seq(0.01,0.99,length = result_all[[1]]$setting$t_total)
  j_vec_eval <- seq(0.01,0.99,length = result_all[[1]]$setting$j_total)
  f1_output <- output_organize(result_all = result_all, start_index = 1, vec_eval = j_vec_eval)
  f2_output <- output_organize(result_all = result_all, start_index = fresult_length+1, vec_eval = t_vec_eval)
  f3_output <- output_organize(result_all = result_all, start_index = 2*fresult_length+1, vec_eval = t_vec_eval)
  max_est_equ <- sapply(result_all, "[[", "est_equ")
  df_all <- sapply(result_all, "[[", "setting")
  #combination of f's
  fcomb_index <- which(names(result_all[[1]]) == "Z1 = 1, at j = 2")
  fcombt_CPdat <- data.frame(t = t_vec_eval)
  for (index in 0:5){
    fcomb_in <- sapply(result_all, "[[", fcomb_index + index)
    fcomb_CP <- apply(fcomb_in, 1, mean)
    fcombt_CPdat[, names(result_all[[1]])[[fcomb_index + index]]] <- fcomb_CP
  }
  fcomb_index <- which(names(result_all[[1]]) == "Z1 = 1, at t = 2")
  fcombj_CPdat <- data.frame(j = j_vec_eval)
  for (index in 0:5){
    fcomb_in <- sapply(result_all, "[[", fcomb_index + index)
    fcomb_CP <- apply(fcomb_in, 1, mean)
    fcombj_CPdat[, names(result_all[[1]])[[fcomb_index + index]]] <- fcomb_CP
  }
  
  alpha0_hat <- sapply(result_all, "[[", "alpha0_hat") 
  alpha0_mse <- mean((alpha0_hat - 0.7)^2)
  alpha0_in <- sapply(result_all, "[[", "alpha0_in") 
  alpha0_cp <- mean(alpha0_in)
  alpha0_result <- c(alpha0_mse, alpha0_cp)
  
  result_all[[i_design]] <- list(f1_output, f2_output, f3_output, max_est_equ, fcombt_CPdat, fcombj_CPdat, df_all, alpha0_result)
}

# Make plots -------------------------------------------

beta_labels <- c(TeX(r'($\hat{\beta}_1$)'),
                 TeX(r'($\hat{\beta}_{21}$)'),
                 TeX(r'($\hat{\beta}_{22}$)'))
x_labels <- c("Time", "Time", "Time")
axis_title_size <- 17
axis_text_size <- 14
plot_list <- list()

for (f_name in 1:3){
  CP_plot <- mse_plot <- knots_plot <- c()
  for (i_design in 1:3){
    result <- result_all[[i_design]][[f_name]]
    cp_dat <- result$cp_dat
    i_total <- result$setting$i_total
    t_total <- result$setting$t_total
    j_total <- result$setting$j_total
    eval_length <- dim(cp_dat)[1]
    
    df_all <- result_all[[i_design]][[7]]
    df_all <- as.data.frame(apply(df_all, 1, as.numeric))
    knots_plot <- rbind(knots_plot, data.frame(knots = df_all[, paste0("df", f_name)] - 3, i_total = i_total))
    
    CP_plot <- rbind(CP_plot, data.frame(cp_dat, i_total = i_total))
    mse_plot <- rbind(mse_plot, data.frame(result$ptconv, i_total = i_total, vec_eval = cp_dat$vec_eval))
  }
  #Plot1: CP
  #colors <- hcl.colors(3, palette = "Teal")
  N_colors <- c("20" = "#88C3C8" , "50" = "#498EA4", "100" = "#2A5676")
  cp1 <- CP_plot %>%
    filter(vec_eval > 0.01001) %>%
    ggplot(aes(x = vec_eval, y = fCP_fhat_varfhat, color = factor(i_total))) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.95, linetype = 2, linewidth = 1) + 
    coord_cartesian(ylim = c(0.7, 1)) +
    theme_bw() +
    theme(legend.key.width =  unit(0.8, 'cm'),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 13), 
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)) +
    ylab("CP") +
    xlab(x_labels[f_name])
  
  # if (f_name == 1) {
  #     cp1 <- cp1 + scale_color_manual(name = "Sample size", values = N_colors)
  # } else {
  #     cp1 <- cp1 + scale_color_manual(guide = "none", values = N_colors)
  # }
  cp1 <- cp1 + scale_color_manual(guide = "none", values = N_colors)
  
  
  #Plot2: MSE
  mse1 <- ggplot(mse_plot, aes(x = vec_eval, y = mse_f_fhat, color = factor(i_total))) +
    geom_line(linewidth = 1) +
    coord_cartesian(ylim = c(0, 0.15)) +
    theme_bw() +
    theme(legend.key.width =  unit(0.8, 'cm'),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 13), 
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)) +
    ylab("MSE") +
    xlab(x_labels[f_name])
  
  mse1 <- mse1 + scale_color_manual(guide = "none", values = N_colors)
  
  
  #Plot3: selected number of knots
  knots_plot_summarized <-  knots_plot %>%
    group_by(i_total) %>%
    table %>%
    as_tibble %>% 
    mutate(knots = as.numeric(knots),
           i_total = as.numeric(i_total))
  
  freq1 <- ggplot(knots_plot_summarized, aes(x = knots, y = n, fill = factor(i_total))) + 
    geom_bar(position="dodge", stat="identity") +
    #geom_freqpoly() +
    theme_bw() +
    theme(legend.key.width =  unit(0.8, 'cm'),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 13), 
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)) +
    ylab("Frequency") +
    xlab("Selected number of knots") +
    coord_cartesian(xlim = c(3,15))
  
  freq1 <- freq1 + scale_fill_manual(guide = "none", values = N_colors)
  
  p1 <- ggarrange(mse1, cp1, freq1, nrow = 1, common.legend = TRUE)
  p1 <- annotate_figure(p1, top = text_grob(beta_labels[f_name], 
                                            face = "bold", size = 15))
  
  plot_list[[f_name]] <- p1
}

ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 1, nrow = 3, legend = "bottom")
ggsave("simulation_result_varyingN_complex_slide.pdf", width = 12, height = 7)

