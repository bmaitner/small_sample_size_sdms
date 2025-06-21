library(disdat)
library(np)
library(AUC)
library(rvinecopulib)
library(pROC)
library(lemon)
library(S4DM)
library(kernlab)
library(tidyverse)
library(sf)
library(DescTools)
library(foreach)
library(doParallel)
source("R/evaluate_disdat.R")

#Select pnp modules to consider (as both numerator and denominator)

  pnp_components <- c("rangebagging",
                      "gaussian",
                      "kde",
                      "vine",
                      "lobagoc")
  
  
  
# Make full set of hybrid models to consider
  models_to_evaluate <- NULL
  
  for(i in pnp_components){
  for(j in pnp_components){
    models_to_evaluate <- 
      rbind(models_to_evaluate,
          data.frame(presence_method = i,
                     background_method = j)) 

  }}

  models_to_evaluate <-
  rbind(models_to_evaluate,
        data.frame(presence_method = pnp_components,
                   background_method = "none"))
  
# Iterate through models
  
full_model_outputs <- NULL
fold_model_outputs <- NULL

# Specify temp file

tempfile_full <- "outputs/temp_bakeoff_output_full.rds"
tempfile_fold <- "outputs/temp_bakeoff_output_fold.rds"
  


  for(i in 1:nrow(models_to_evaluate)){
    
    # If the temporary output files exist, check to see what they contain
    
      if(file.exists(file.path(tempfile_fold))){
        fold_model_outputs <- readRDS(file.path(tempfile_fold))
      }
    
      if(file.exists(file.path(tempfile_full))){
        full_model_outputs <- readRDS(file.path(tempfile_full))
      }
    
    # If model has already been done, move on to next
    
      presence_method_i <- models_to_evaluate$presence_method[i]
      background_method_i <- models_to_evaluate$background_method[i]
      
      message("presence method: ", presence_method_i,
              "; background method: ",background_method_i)    
  
      if( any(fold_model_outputs$pres_method == presence_method_i &
              fold_model_outputs$bg_method == background_method_i) &
          any(full_model_outputs$pres_method == presence_method_i &
              full_model_outputs$bg_method == background_method_i)){next}
      
    set.seed(2005) # The year Transformers: The Movie is set.
    
    model_i <- 
    evaluate_disdat(presence_method = models_to_evaluate$presence_method[i],
                    background_method = models_to_evaluate$background_method[i],
                    verbose = TRUE,
                    ratio_method = NULL,
                    quantile = 0.05,
                    ncl = 5,
                    record_predictions = TRUE,
                    predictions_folder = "outputs/model_predictions/")
    
    full_model_outputs <- rbind(full_model_outputs,
                                data.frame(pres_method = models_to_evaluate$presence_method[i],
                                           bg_method = models_to_evaluate$background_method[i],
                                           model_i$full_model_stats))
    
    fold_model_outputs <- rbind(fold_model_outputs,
                                data.frame(pres_method = models_to_evaluate$presence_method[i],
                                           bg_method = models_to_evaluate$background_method[i],
                                           model_i$fold_model_stats))
    
    # Save temporary files

      fold_model_outputs %>%
        saveRDS(file = file.path(tempfile_fold))
    
      full_model_outputs %>%
        saveRDS(file.path(tempfile_full))
  }
  
# save outputs as an RDS object (since it takes so long to re-run)

# saveRDS(object = full_model_outputs,
#         file = "outputs/bake_off_pnp_full_model_outputs.RDS")
# 
# saveRDS(object = fold_model_outputs,
#         file = "outputs/bake_off_pnp_fold_model_outputs.RDS")

fold_model_outputs <- readRDS("outputs/bake_off_pnp_fold_model_outputs.RDS")
full_model_outputs <- readRDS("outputs/bake_off_pnp_full_model_outputs.RDS")

# ggplots of model stats (facet grid of numerator and denominator)  
  
  
library(ggplot2)  
##################################

dr_models_to_evaluate <- c("ulsif","rulsif","maxnet")

# Iterate through models

# Iterate through models

full_model_outputs_dr <- NULL
fold_model_outputs_dr <- NULL

# Specify temp file

tempfile_full_dr <- "outputs/temp_bakeoff_output_full_dr.rds"
tempfile_fold_dr <- "outputs/temp_bakeoff_output_fold_dr.rds"

for(i in 1:length(dr_models_to_evaluate)){
  
  # If the temporary output files exist, check to see what they contain
  
  if(file.exists(file.path(tempfile_fold_dr))){
    fold_model_outputs_dr <- readRDS(file.path(tempfile_fold_dr))
  }
  
  if(file.exists(file.path(tempfile_full_dr))){
    full_model_outputs_dr <- readRDS(file.path(tempfile_full_dr))
  }
  

  message("density ratio method: ", dr_models_to_evaluate[i])    
  
  if( any(fold_model_outputs_dr$ratio_method == dr_models_to_evaluate[i]) ){next}
  
  set.seed(2005) # The year Transformers: The Movie is set.

  model_i <- 
    evaluate_disdat(ratio_method = dr_models_to_evaluate[i])
  
  
  full_model_outputs_dr <- rbind(full_model_outputs_dr,
                              data.frame(ratio_method = dr_models_to_evaluate[i],
                                         model_i$full_model_stats))
  
  fold_model_outputs_dr <- rbind(fold_model_outputs_dr,
                              data.frame(ratio_method = dr_models_to_evaluate[i],
                                         model_i$fold_model_stats))
  
  
  # Save temporary files
  
  fold_model_outputs_dr %>%
    saveRDS(file = file.path(tempfile_fold_dr))
  
  full_model_outputs_dr %>%
    saveRDS(file.path(tempfile_full_dr))
  
}


# save outputs as an RDS object (since it takes so long to re-run)

# saveRDS(object = full_model_outputs_dr,
#         file = "outputs/bake_off_dr_full_model_outputs.RDS")
# 
# saveRDS(object = fold_model_outputs_dr,
#         file = "outputs/bake_off_dr_fold_model_outputs.RDS")


fold_model_outputs_dr <- readRDS("outputs/bake_off_dr_fold_model_outputs.RDS")
full_model_outputs_dr <- readRDS("outputs/bake_off_dr_full_model_outputs.RDS")

##################################

#Combing dr and pnp results

#Make sure colnames match up

fold_model_outputs %>% 
  mutate(method = paste(pres_method,"/",bg_method),
         ratio_method = NA) -> fold_model_outputs

full_model_outputs %>% 
  mutate(method = paste(pres_method,"/",bg_method),
         ratio_method = NA) -> full_model_outputs


fold_model_outputs_dr %>% 
  mutate(method = ratio_method,
         pres_method = NA,
         bg_method = NA) -> fold_model_outputs_dr

full_model_outputs_dr %>% 
  mutate(method = ratio_method,
         pres_method = NA,
         bg_method = NA) -> full_model_outputs_dr


#Re-arrange columns and merge

full_model_outputs_dr <- full_model_outputs_dr[colnames(full_model_outputs)]
fold_model_outputs_dr <- fold_model_outputs_dr[colnames(fold_model_outputs)]


fold_model_output_all <- rbind(fold_model_outputs,fold_model_outputs_dr)
full_model_output_all <- rbind(full_model_outputs,full_model_outputs_dr)


#saveRDS(object = full_model_output_all,file = "outputs/full_model_output_all.RDS")
#saveRDS(object = fold_model_output_all,file = "outputs/fold_model_output_all.RDS")


########################

library(tidyverse)

library(ggrepel)
library(grid)


temp<-
full_model_output_all %>%
  #filter(n_presence <= 20)%>%
  filter(method %in% c("kde / none","vine / none","gaussian / none",
                      "maxnet","ulsif","rulsif",
                      "rangebagging / none","lobagoc / none"))%>%
  mutate(method = gsub(pattern = " / none",replacement = "",x=method))%>%
  mutate(method = gsub(pattern = "maxnet",replacement = "Maxnet",x=method))%>%
  mutate(method = gsub(pattern = "vine",replacement = "Vine",x=method))%>%
  mutate(method = gsub(pattern = "gaussian",replacement = "Gaussian",x=method))%>%
  mutate(method = gsub(pattern = "kde",replacement = "KDE",x=method))%>%
  mutate(method = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=method))%>%
  mutate(method = gsub(pattern = "ulsif",replacement = "uLSIF",x=method))%>%
  mutate(method = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=method))%>%
  mutate(sens_spec_ratio = pa_specificity-pa_sensitivity) %>%
  group_by(method)%>%
  summarise(sens_spec_ratio = median(na.omit(sens_spec_ratio)))%>%
  ggplot(mapping = aes(x=sens_spec_ratio,y=0,label = method))+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_text_repel(max.overlaps = 20,box.padding = 1)+
  xlab("Sensitivity - Specificity")+
  ylab(NULL)+
  ylim(c(-.5,.5))+
  xlim(c(-1,1))+
  theme_bw()+
  guides(color="none")+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        line = element_blank(),
        axis.title =element_text(size = 15))
  
temp

text_high_spec <- textGrob("High Specificity\nPresences Correct\nAssumes good sampling",)
text_high_sens <- textGrob("High Sensitivity\nAbsences Correct\nAssumes poor sampling")


sens_spec_gradient<-
temp+
  annotation_custom(text_high_sens,xmin=-1,xmax=-0.5,ymin=-1,ymax=-.65)+
  annotation_custom(text_high_spec,xmin=0.5,xmax=1,ymin=-1,ymax=-.65)+
  theme(plot.margin = unit(c(1,1,1,1), "lines")) + #top,right,bottom,left
  theme(plot.margin = unit(c(1,1,6,1), "lines")) + #top,right,bottom,left
  coord_cartesian(ylim=c(-.4,.4), clip="off")

sens_spec_gradient

ggsave(plot = sens_spec_gradient,filename = "figures/sens_spec_gradient.jpg",
       width = 6,height = 4,units = "in",dpi = 600)

ggsave(plot = sens_spec_gradient,filename = "figures/sens_spec_gradient.svg",
       width = 6,height = 4,units = "in",dpi = 600)

################

#How many time does lobagoc beat maxnet?

full %>%
  filter(method %in% c("rangebagging / none","maxnet")) %>%
  group_by(species)%>%
  arrange(species,pa_AUC)%>%
  slice_head(n = 1)%>%
  ungroup()%>%
  select(method)%>%
  table()

#88/(88+138) #39%


full %>%
  filter(method %in% c("lobagoc / none","maxnet")) %>%
  group_by(species)%>%
  arrange(species,desc(pa_AUC))%>%
  select(species,method,pa_AUC)%>%
  slice_head(n = 1)%>%
  ungroup()%>%
  select(method)%>%
  table()

  #45/(45+181) #20%

##########################################################################


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(grid)


full_model_output_all <- readRDS(file = "outputs/full_model_output_all.RDS")

full_model_output_all %>%
  #filter(n_presence <= 20)%>%
  filter(method %in% c("kde / none","vine / none","gaussian / none",
                       "maxnet","ulsif","rulsif",
                       "rangebagging / none","lobagoc / none"))%>%
  mutate(method = gsub(pattern = " / none",replacement = "",x=method))%>%
  mutate(method = gsub(pattern = "maxnet",replacement = "Maxnet",x=method))%>%
  mutate(method = gsub(pattern = "vine",replacement = "Vine",x=method))%>%
  mutate(method = gsub(pattern = "gaussian",replacement = "Gaussian",x=method))%>%
  mutate(method = gsub(pattern = "kde",replacement = "KDE",x=method))%>%
  mutate(method = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=method))%>%
  mutate(method = gsub(pattern = "ulsif",replacement = "uLSIF",x=method))%>%
  mutate(method = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=method)) %>%  
  #mutate(sens_spec_ratio = pa_specificity-pa_sensitivity) %>%
  group_by(method)%>%
  #summarise(sens_spec_ratio = median(na.omit(sens_spec_ratio)))%>%
  summarise(mean_sensitivity = mean(na.omit(pa_sensitivity)),
            mean_specificity = mean(na.omit(pa_specificity)),
            ci_low_sensitivity = quantile(x = pa_sensitivity,
                                          probs = 0.25,
                                          na.rm=TRUE),
            ci_high_sensitivity = quantile(x = pa_sensitivity,
                                           probs = 0.75,
                                           na.rm=TRUE),
            ci_low_specificity = quantile(x = pa_specificity,
                                          probs = 0.25,
                                          na.rm=TRUE),
            ci_high_specificity = quantile(x = pa_specificity,
                                           probs = 0.75,
                                           na.rm=TRUE)
            
  ) %>%
  ggplot(mapping = aes(x=mean_sensitivity,
                       y=mean_specificity,
                       label = method,
                       color = method))+
  geom_errorbar(mapping = aes(ymin=ci_low_specificity,
                              ymax=ci_high_specificity),
                linewidth = 1)+
  geom_errorbarh(mapping = aes(xmin=ci_low_sensitivity,
                               xmax=ci_high_sensitivity),
                 linewidth = 1)+
  geom_point(size=3)+
  scale_x_continuous(expand = c(0, 0),limits = c(0,1.2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.2))+
  ylab("Specificity")+
  xlab("Sensitivity")+
  theme_bw()+
  guides(color="none")+ 
  geom_text_repel(max.overlaps = 20,
                  box.padding = 3,
                  point.padding = 0)->svs

ggsave(plot = svs,
       filename = "figures/sensitivity_v_specificity_pres_only.jpg",
       width = 5,
       height = 5,
       units = "in",dpi = 600)

ggsave(plot = svs,
       filename = "figures/sensitivity_v_specificity_pres_only.svg",
       width = 5,
       height = 5,
       units = "in",dpi = 600)


#######################################################################

# the subtypes together model all 

full_model_output_all %>%
  mutate(method = case_when(is.na(pres_method) ~ ratio_method,
                            !is.na(pres_method) ~pres_method)) %>%
  mutate(method = gsub(pattern = " / none",replacement = "",x=method))%>%
  mutate(method = gsub(pattern = "maxnet",replacement = "Maxnet",x=method))%>%
  mutate(method = gsub(pattern = "vine",replacement = "Vine",x=method))%>%
  mutate(method = gsub(pattern = "gaussian",replacement = "Gaussian",x=method))%>%
  mutate(method = gsub(pattern = "kde",replacement = "KDE",x=method))%>%
  mutate(method = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=method))%>%
  mutate(method = gsub(pattern = "ulsif",replacement = "uLSIF",x=method))%>%
  mutate(method = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=method))%>%
  group_by(method)%>%
  #summarise(sens_spec_ratio = median(na.omit(sens_spec_ratio)))%>%
  summarise(mean_sensitivity = mean(na.omit(pa_sensitivity)),
            mean_specificity = mean(na.omit(pa_specificity)),
            ci_low_sensitivity = quantile(x = pa_sensitivity,
                                          probs = 0.25,
                                          na.rm=TRUE),
            ci_high_sensitivity = quantile(x = pa_sensitivity,
                                           probs = 0.75,
                                           na.rm=TRUE),
            ci_low_specificity = quantile(x = pa_specificity,
                                          probs = 0.25,
                                          na.rm=TRUE),
            ci_high_specificity = quantile(x = pa_specificity,
                                           probs = 0.75,
                                           na.rm=TRUE)
            
  ) %>%
  ggplot(mapping = aes(x=mean_sensitivity,
                       y=mean_specificity,
                       label = method,
                       color = method))+
  geom_errorbar(mapping = aes(ymin=ci_low_specificity,
                              ymax=ci_high_specificity),
                linewidth = 1)+
  geom_errorbarh(mapping = aes(xmin=ci_low_sensitivity,
                               xmax=ci_high_sensitivity),
                 linewidth = 1)+
  geom_point(size=3)+
  scale_x_continuous(expand = c(0, 0),limits = c(0,1.2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,1.2))+
  ylab("Specificity")+
  xlab("Sensitivity")+
  theme_bw()+
  guides(color="none")+ 
  geom_text_repel(max.overlaps = 20,
                  box.padding = 3,
                  point.padding = 0)->svs2

plot(svs2)  

ggsave(plot = svs2,
       filename = "figures/sensitivity_v_specificity.jpg",
       width = 5,
       height = 5,
       units = "in",dpi = 600)

ggsave(plot = svs2,
       filename = "figures/sensitivity_v_specificity.svg",
       width = 5,
       height = 5,
       units = "in",dpi = 600)


