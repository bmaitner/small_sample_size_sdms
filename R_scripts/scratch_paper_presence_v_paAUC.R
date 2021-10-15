library(tidyverse)
library(ggplot2)
library(gganimate)

temp_rarified <- readRDS("outputs/temp_rarified.RDS")

temp_rarified %>%
  #group_by(model) %>%
  filter(model %in% c("maxnet","gaussian/kde","ulsif","kde/kde","rangebagging/none","gaussian/gaussian","kde/none"))%>%
  #filter(n_presence < 50)%>%
  #ggplot( mapping = aes(y= pa_AUC, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_sensitivity, x = n_presence^.5, color = model))+
  ggplot( mapping = aes(y= pa_precision, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_specificity, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_prediction_accuracy, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_pAUC_sensitivity, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_pAUC_specificity, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_kappa, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_correlation, x = n_presence^.5, color = model))+
  geom_point(alpha=1)+geom_smooth(se = T)+
  #ylab("presence-absence AUC")+
  xlab("Number of Presences ^0.5")

temp_rarified %>%
  #group_by(model) %>%
  filter(model %in% c("maxnet","gaussian/kde","ulsif","kde/kde","rangebagging/none","gaussian/gaussian"))%>%
  #filter(n_presence < 50)%>%
  ggplot( mapping = aes(y= runtime, x = n_presence^0.5, color = model))+
  #ggplot( mapping = aes(y= pa_sensitivity, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_prediction_accuracy, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_pAUC_sensitivity, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_pAUC_specificity, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_kappa, x = n_presence^.5, color = model))+
  #ggplot( mapping = aes(y= pa_correlation, x = n_presence^.5, color = model))+
  geom_point(alpha=1)+geom_smooth(se = T)+
  #ylab("presence-absence AUC")+
  xlab("Number of Presences ^0.5")




?gganimate
?gg_animate

anim <- ggplot(mtcars, aes(mpg, disp)) +
  transition_states(gear, transition_length = 2, state_length = 1) +
  enter_fade() +
  exit_fade()

animate(anim)
animate(anim, fps = 20, duration = 15)



