rm(list = ls())

library(tidyverse)
library(xtable)
library(scico)

jbslot60 <- readRDS("jbslot_public_60min.RDS")
total_T <- 60

## one heatmap for each user ##

# Reconstruct user-specific data frame, with decision points ordered first by send/not send (optional: then by clustering result).
# This way, the heatmap may be easier for us to see patterns.

## heatmap without clustering

##categorical treatment using log(steps)
jbslot60 %>% 
  filter(user == 1) %>% 
  group_by(user) %>%
  do({
    tmp_active <- matrix(.$steps.log[.$send.active == 1], nrow = 60, byrow = FALSE)
    tmp_seden <- matrix(.$steps.log[.$send.sedentary == 1], nrow = 60, byrow = FALSE)
    tmp0 <- matrix(.$steps.log[.$send == 0], nrow = 60, byrow = FALSE)
    
    n_dp <- nrow(.) / total_T
    n_send0 <- sum(.$send == 0)
    n_send_active <- sum(.$send.active == 1)
    n_send_seden <- sum(.$send.sedentary == 1)
    print(n_dp)
    print(total_T)
    print(length(tmp0))
    print(length(tmp_seden))
    print(length(tmp_active))
      
    user_clus <- data.frame(
      dp_pseudo = rep(0:(n_dp - 1), each = total_T),
      min.after.decision = rep(1:total_T, n_dp), 
      steps.log = c(as.vector(tmp0), as.vector(tmp_seden), as.vector(tmp_active)))
    #separate them
    user_clus$dp_pseudo[(n_send0 + 1):(n_send0 + n_send_seden)] <- user_clus$dp_pseudo[(n_send0 + 1):(n_send0 + n_send_seden)] + 4
    user_clus$dp_pseudo[(n_send0 + n_send_seden + 1):nrow(user_clus)] <- 
      user_clus$dp_pseudo[(n_send0 + n_send_seden + 1):nrow(user_clus)] + 8
    
    print(head(user_clus))
    
    p <- ggplot(data = user_clus, mapping = aes(x = min.after.decision,
                                                y = dp_pseudo,
                                                fill = steps.log)) + 
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=20),
            legend.text = element_text(size=18),
            legend.title = element_text(size=20),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      # scale_fill_continuous(name = "log(steps)", type = "viridis") +
      scale_fill_scico(name = "log(steps)", palette = 'acton', direction = -1) +
      xlab("minute") +
      geom_tile() + ylab("")
    
    ggsave(p, filename = paste0("user", .$user[1], "_heatmap_nogray_categorical_logsteps.png"), width = 11.1, height = 9.78) 
    
    invisible(.) # this line seems important; otherwise all kinds of error
  })


##categorical treatment using steps
jbslot60 %>% 
  filter(user == 1) %>% 
  group_by(user) %>%
  do({
    tmp_active <- matrix(.$steps[.$send.active == 1], nrow = 60, byrow = FALSE)
    tmp_seden <- matrix(.$steps[.$send.sedentary == 1], nrow = 60, byrow = FALSE)
    tmp0 <- matrix(.$steps[.$send == 0], nrow = 60, byrow = FALSE)
    
    n_dp <- nrow(.) / total_T
    n_send0 <- sum(.$send == 0)
    n_send_active <- sum(.$send.active == 1)
    n_send_seden <- sum(.$send.sedentary == 1)
    
    user_clus <- data.frame(
      dp_pseudo = rep(0:(n_dp - 1), each = total_T),
      min.after.decision = rep(1:total_T, n_dp), 
      steps.log = c(as.vector(tmp0), as.vector(tmp_seden), as.vector(tmp_active)))
    #separate them
    user_clus$dp_pseudo[(n_send0 + 1):(n_send0 + n_send_seden)] <- user_clus$dp_pseudo[(n_send0 + 1):(n_send0 + n_send_seden)] + 4
    user_clus$dp_pseudo[(n_send0 + n_send_seden + 1):nrow(user_clus)] <- 
      user_clus$dp_pseudo[(n_send0 + n_send_seden + 1):nrow(user_clus)] + 8
    
    print(head(user_clus))
    
    p <- ggplot(data = user_clus, mapping = aes(x = min.after.decision,
                                                y = dp_pseudo,
                                                fill = steps.log)) + 
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=18),
            axis.title=element_text(size=20),
            legend.text = element_text(size=18),
            legend.title = element_text(size=20),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      # scale_fill_continuous(name = "log(steps)", type = "viridis") +
      scale_fill_scico(name = "steps", palette = 'acton', direction = -1) +
      xlab("minute") +
      geom_tile() + ylab("")
    
    ggsave(p, filename = paste0("user", .$user[1], "_heatmap_nogray_categorical_steps.png"), width = 11.1, height = 9.78) 
    
    invisible(.) # this line seems important; otherwise all kinds of error
  })

# user's step count plot --------------------------------------------------
## making plot for individual's step counts 
title_size <- 15

p_2 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 1) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
#  ggtitle("Walking suggestion \n(at decision point 2)") +
  ggtitle("DP 2 (Walking suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))

p_12 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 11) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
  ggtitle("DP 12 (Anti-sedentary suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))

p_31 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 30) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
#  ggtitle("Walking suggestion \n(at decision point 31)") +
  ggtitle("DP 31 (Walking suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))

p_61 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 60) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
#  ggtitle("Walking suggestion \n(at decision point 61)") +
  ggtitle("DP 61 (Walking suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))

p_104 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 103) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
#  ggtitle("No suggestion \n(at decision point 104)") +
  ggtitle("DP 104 (No suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))

p_145 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 144) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
#  ggtitle("Anti-sedentary suggestion \n(at decision point 145)") +
  ggtitle("DP 145 (Anti-sedentary suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))


p_173 <- jbslot60 %>% 
  filter(user == 1, decision.index.nogap == 172) %>% 
  ggplot(aes(x = min.after.decision, y = steps)) + geom_point(size = 0.5) + geom_line() +
  coord_cartesian(ylim = c(0, 200)) +
  xlab("minute") +
  ylab("step count") +
#  ggtitle("No suggestion \n(at decision point 173)") +
  ggtitle("DP 173 (No suggestion)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=title_size))

library(ggpubr)
p <- ggarrange(p_2, p_31, 
               p_12, p_145, 
               p_104, p_173, nrow = 3, ncol = 2)
ggsave(p, filename = "user1 step count plot.png", width = 8, height = 8.5)



