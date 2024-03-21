library(disdat)
library(np)
library(AUC)
library(rvinecopulib)
library(pROC)
library(lemon)
library(pbsdm)
library(tidyverse)
library(sf)
library(DescTools)
library(foreach)
library(doParallel)
source("R/evaluate_disdat.R")

#Select pnp modules to consider (as both numerator and denominator)
  pnp_components <- c("rangebagging",
                      "gaussian",
                      "kde")
  
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
                    quantile = 0.05)
    
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

##################################
#Plot: pres x bg method PA AUC histograms


full_model_outputs %>%
  group_by(pres_method,bg_method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(presence_category = 
           case_when(n_presence < min_threshold ~ "few",
                     n_presence > max_threshold ~ "many",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "intermediate"
                       ))%>%ungroup%>%
  mutate(bg_method = fct_relevel(bg_method,"none","rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"),
         presence_category = fct_relevel(presence_category, "few","intermediate","many")) %>%
  #filter(point_category != "typical")%>%
  ggplot( mapping = aes(x= pa_AUC, fill = presence_category))+
  #ggplot( mapping = aes(x= pa_AUC, fill = pres_method))+
  #geom_histogram(position = "identity", alpha = 0.5)+
  geom_density(alpha=0.5)+
  geom_vline(data = . %>%
      group_by(bg_method,pres_method) %>%
      summarise(line = median(pa_AUC,na.rm=T)),
    mapping = aes(xintercept = line)
  )+
  facet_rep_grid(pres_method ~ bg_method,repeat.tick.labels = T,
             labeller = labeller(
               pres_method = function(x){paste("Presence: \n", x)},
               bg_method = function(x){paste("Background: \n", x)})
             )+

  geom_vline(
    data = . %>%
      group_by(bg_method, pres_method) %>%
      summarise(line = quantile(pa_AUC,
                                probs = .75,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  geom_vline(
    data = . %>%
      group_by(bg_method, pres_method) %>%
      summarise(line = quantile(pa_AUC,
                                probs = .25,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  xlab("Presence-Absence AUC")+
  labs(fill = "Presence category")+
  ylab("Density")


full_model_outputs_dr %>%
  group_by(ratio_method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(presence_category = 
           case_when(n_presence < min_threshold ~ "few",
                     n_presence > max_threshold ~ "many",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "intermediate"
           ))%>%ungroup%>%
  mutate(ratio_method = fct_relevel(ratio_method,"ulsif","rulsif","maxnet"),
         presence_category = fct_relevel(presence_category, "few","intermediate","many")) %>%
  #filter(point_category != "typical")%>%
  ggplot( mapping = aes(x= pa_AUC, fill = presence_category))+
  #ggplot( mapping = aes(x= pa_AUC, fill = pres_method))+
  #geom_histogram(position = "identity", alpha = 0.5)+
  geom_density(alpha=0.5)+
  geom_vline(data = . %>%
               group_by(ratio_method) %>%
               summarise(line = median(pa_AUC, na.rm = T)),
             mapping = aes(xintercept = line)
  )+
  facet_rep_wrap("ratio_method",repeat.tick.labels = T)+
  geom_vline(
    data = . %>%
      group_by(ratio_method) %>%
      summarise(line = quantile(pa_AUC,
                                probs = .75,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  geom_vline(
    data = . %>%
      group_by(ratio_method) %>%
      summarise(line = quantile(pa_AUC,
                                probs = .25,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  xlab("Presence-Absence AUC")+
  labs(fill = "Presence category")+
  ylab("Density")+xlim(c(0.4,1))












##############################################  
  
#Figure 3 in Wisz et al.
  

#median auc vs sample size, color = method
  
colnames(full_model_output_all)
  
full_model_output_all %>%
  group_by(species) %>%
  mutate(mn_rel_pa_AUC = pa_AUC-.data$pa_AUC[[which(.data$method == "maxnet")]])%>%
  #filter(method %in% c("gaussian / gaussian", "maxnet", "gaussian / kde")) %>%  
ggplot(mapping = aes(x = log10(n_presence),
                     y = mn_rel_pa_AUC,
                     color = method)) +
  geom_point()+geom_smooth()
  
full_model_output_all %>%
  group_by(species) %>%
  mutate(mn_rel_pa_AUC = pa_AUC-.data$pa_AUC[[which(.data$method=="maxnet")]])%>%
  filter(bg_method %in% c("rangebagging") | pres_method %in% c("rangebagging") | method %in% c("maxnet")) %>%  
  filter(n_presence < 50) %>%
  ggplot(mapping = aes(x = log10(n_presence),
                       y =pa_AUC,
                       color = method)) +
  #geom_point()+
  geom_smooth(se = F)



full_model_output_all %>%
  group_by(species) %>%
  mutate(mn_rel_pa_AUC = pa_AUC-.data$pa_AUC[[which(.data$method=="maxnet")]])%>%
  filter(bg_method %in% c("gaussian") | pres_method %in% c("gaussian") | method %in% c("maxnet")) %>%  
  filter(n_presence < 50) %>%
    ggplot(mapping = aes(x = log10(n_presence),
                       y =pa_AUC,
                       color = method)) +
  #geom_point()+
  geom_smooth(se = F)

full_model_output_all %>%
  group_by(species) %>%
  mutate(mn_rel_pa_AUC = pa_AUC-.data$pa_AUC[[which(.data$method=="maxnet")]])%>%
  filter(bg_method %in% c("kde") | pres_method %in% c("kde") | method %in% c("maxnet")) %>%  
  filter(n_presence < 100) %>%
  ggplot(mapping = aes(x = log10(n_presence),
                       y =pa_AUC,
                       color = method)) +
  #geom_point()+
  geom_smooth(se = F)





##############################################

full_model_outputs %>%
  group_by(pres_method,bg_method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%
  mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde")) %>%
  filter(point_category != "typical")%>%
  ggplot( mapping = aes(x= pa_AUC, fill = point_category))+
  #geom_histogram(position = "identity", alpha = 0.5)+
  geom_density(alpha=0.5)+
  geom_vline(data = . %>%
               group_by(bg_method,pres_method) %>%
               summarise(line = median(pa_AUC,na.rm=T)),
             mapping = aes(xintercept = line)
  )+
  facet_grid(pres_method ~ bg_method,
             labeller = labeller(
               pres_method = function(x){paste("Presence: \n", x)},
               bg_method = function(x){paste("Background: \n", x)})
  )+
  
  geom_vline(
    data = . %>%
      group_by(bg_method, pres_method) %>%
      summarise(line = quantile(pa_AUC,
                                probs = .75,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  geom_vline(
    data = . %>%
      group_by(bg_method, pres_method) %>%
      summarise(line = quantile(pa_AUC,
                                probs = .25,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  xlab("presence-absence AUC")

+
  theme(legend.position = "none")


###############################################

full_model_outputs %>%
  mutate(bg_method = fct_relevel(bg_method,"none","rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"none","rangebagging","gaussian","kde")) %>%
  # mutate(bg_method = paste("Background: ", .data$bg_method, sep = ""),
  #        pres_method = paste("Presence: ", .data$pres_method, sep = "")) %>%
  ggplot( mapping = aes(x= full_AUC, fill = pres_method))+
  geom_histogram(position = "identity", alpha = 0.5)+
  geom_vline(data = . %>%
               group_by(bg_method,pres_method) %>%
               summarise(line = median(full_AUC,na.rm=T)),
             mapping = aes(xintercept = line)
  )+
  facet_grid(pres_method ~ bg_method,
             labeller = labeller(
               pres_method = function(x){paste("Presence: \n", x)},
               bg_method = function(x){paste("Background: \n", x)})
  )+
  
  geom_vline(
    data = . %>%
      group_by(bg_method, pres_method) %>%
      summarise(line = quantile(full_AUC,
                                probs = .75,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  geom_vline(
    data = . %>%
      group_by(bg_method, pres_method) %>%
      summarise(line = quantile(full_AUC,
                                probs = .25,na.rm=T)),
    mapping = aes(xintercept = line),
    linetype=2
  )+
  xlab("AUC")+
  theme(legend.position = "none")






library(tidyverse)
full_model_outputs %>%
  group_by(bg_method,pres_method) %>%
  summarise(mean_auc = mean(na.omit(full_AUC)),
            mean_pa_auc = mean(na.omit(pa_AUC)),
            mean_pa_pAUC_sensitivity = mean(na.omit(pa_pAUC_sensitivity)),
            mean_pa_pAUC_specificity = mean(na.omit(pa_pAUC_specificity)),
            mean_pa_prediction_accuracy = mean(na.omit(pa_prediction_accuracy)),
            mean_pa_sensitivity = mean(na.omit(pa_sensitivity)),
            mean_pa_specificity = mean(na.omit(pa_specificity)),
            mean_pa_pAUC_specificity = mean(na.omit(pa_specificity)),
            mean_pa_correlation = mean(na.omit(pa_correlation)),
            mean_pa_kappa = mean(na.omit(pa_kappa)),
            median_pa_auc =median(na.omit(pa_AUC))
            )%>%
  write.csv(file = "outputs/pa_scores.csv")


full_model_output_all %>%
  group_by(method) %>%
  summarise(mean_auc = mean(na.omit(full_AUC)),
            mean_pa_auc = mean(na.omit(pa_AUC)),
            mean_pa_pAUC_sensitivity = mean(na.omit(pa_pAUC_sensitivity)),
            mean_pa_pAUC_specificity = mean(na.omit(pa_pAUC_specificity)),
            mean_pa_prediction_accuracy = mean(na.omit(pa_prediction_accuracy)),
            mean_pa_sensitivity = mean(na.omit(pa_sensitivity)),
            mean_pa_specificity = mean(na.omit(pa_specificity)),
            mean_pa_correlation = mean(na.omit(pa_correlation)),
            mean_pa_kappa = mean(na.omit(pa_kappa)),
            median_pa_auc =median(na.omit(pa_AUC))
  )%>%
  write.csv(file = "outputs/pa_scores_dr_and_pnp.csv")



colnames(full_model_outputs)

# full_model_outputs$pres_method <- factor(full_model_outputs$pres_method,levels = c("rangebagging","gaussian","kde"),ordered = T)
# full_model_outputs$bg_method <- factor(full_model_outputs$bg_method,levels = c("rangebagging","gaussian","kde"),ordered = T)
# full_model_outputs$pres_method <- factor(full_model_outputs$pres_method,ordered = F)
# full_model_outputs$bg_method <- factor(full_model_outputs$bg_method,ordered = F)


#####################################################################################

# Presence-background model performance

library(lme4)
library(betareg)
library(bbmle)
library(glmmTMB)

full_model_outputs %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  group_by(pres_method,bg_method) %>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%ungroup%>%
  mutate(bg_method = fct_relevel(bg_method,"none","rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"),
         point_category= fct_relevel(point_category, "low","typical","high"))%>%
  ungroup() %>%
  glmmTMB(formula = pa_AUC ~ pres_method + bg_method + pres_method*bg_method +
                point_category + pres_method*point_category
              +(1|species),
              family = list(family = "beta", link = "logit")) -> pb_betareg_out 

#3-way pres x bg x point category not supported

summary(pb_betareg_out)
summary(pb_betareg_out)$coefficients$cond %>% write.csv(file = "outputs/pb_beta_model_summary.csv",row.names = T)

aov(pb_betareg_out)

#difflsmeans(pb_betareg_out,test.effs="Group")
#difflsmeans(pnp_lmer, test.effs = "Group", ddf="Kenward-Roger")

# contrasts(pnp_lmer)
# contrasts(pb_betareg_out)

# glht(model = pnp_lmer, linfct = mcp(group_var = "pres_method"))
# glht(model = pb_betareg_out, linfct = mcp(group_var = "pres_method"))

library(emmeans)
?emmeans

emmeans::emmeans(object = pb_betareg_out,specs = "pres_method")
emmeans::emmeans(object = pb_betareg_out,specs = point_category ~ pres_method * bg_method )


glht(emmeans::emmeans(object = pb_betareg_out,specs = "pres_method"))

emmeans::emmip(object = pb_betareg_out, pres_method ~ bg_method,CIs=TRUE)
emmeans::emmip(object = pb_betareg_out, pres_method * bg_method ~ point_category ,CIs=TRUE,type="response")

emmeans::emmip_ggplot(emms = emmeans::emmip(object = pb_betareg_out, pres_method * bg_method ~ point_category ,CIs=TRUE,type="response",plotit=FALSE),
                      xlab = "Point Category",lty= c(1,2,3,4))

emmeans::emmip_ggplot(emms = emmeans::emmip(object = pb_betareg_out, pres_method * bg_method ~ point_category ,type="response",plotit=FALSE),
                      xlab = "Presence Point Category",
                      linearg = list(linetype=sort(rep(1:12,3)),size=2),
                      tlab = "Presence method\nBackground Method")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



library(multcomp)
library(lmerTest)

##############################################################

#overall model comparison

full_model_output_all %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  group_by(method) %>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%ungroup%>%
  mutate(point_category= fct_relevel(point_category, "low","typical","high"))%>%
  ungroup() %>%
  glmmTMB(formula = pa_AUC ~ method + point_category + method*point_category
          +(1|species),
          family = list(family = "beta", link = "logit")) -> overall_betareg_out

  emmeans::emmip_ggplot(emms = emmeans::emmip(object = overall_betareg_out, method ~ point_category ,type="response",plotit=FALSE),
                        xlab = "Presence Point Category",
                        size=2,
                        tlab = "Method")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  full_model_output_all %>%
    mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
           max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
    group_by(method) %>%
    mutate(point_category = 
             case_when(n_presence < min_threshold ~ "low",
                       n_presence > max_threshold ~ "high",
                       n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
             ))%>%ungroup%>%
    mutate(point_category= fct_relevel(point_category, "low","typical","high"))%>%
    ungroup() %>% filter(method %in% c("maxnet", "ulsif", "rulsif", "gaussian / gaussian", "gaussian / kde", "rangebagging / none",
                                       "kde / kde")) %>%
    glmmTMB(formula = pa_AUC ~ method + point_category + method*point_category
            +(1|species),
            family = list(family = "beta", link = "logit"))-> test
    
    
  emmeans::emmip_ggplot(emms = emmeans::emmip(object=test, method ~ point_category ,type="response",plotit=FALSE,CI=TRUE),
                        xlab = "Presence Point Category",
                        size=2,
                        tlab = "Method")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  

###############################################



full_model_outputs %>%
  group_by(pres_method,bg_method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%
  mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"))%>%
  ungroup()%>%
  group_by(pres_method, bg_method, point_category)%>%
  summarise(mean_pa_AUC = mean(na.omit(pa_AUC)),
            median_pa_AUC = median(na.omit(pa_AUC)))%>%
  #filter(bg_method == "kde")%>%
  filter(!is.na(point_category))%>%
  write.csv(file = "outputs/point_class_pa_scores.csv")


#############################################################

#Statistical models

colnames(full_model_outputs)

# This gives a very different set of results, most likely due to overfitting
# lm_3way <- lm(data = full_model_outputs,formula = full_AUC ~ pres_method + bg_method + pres_method*bg_method + n_presence + pres_method*bg_method*n_presence )
# lm_2way <- lm(data = full_model_outputs,formula = full_AUC ~ pres_method + bg_method + pres_method*bg_method + n_presence )


lm_3way <- lm(data = full_model_outputs,
              formula = pa_AUC ~ pres_method + bg_method + pres_method*bg_method + n_presence + pres_method*bg_method*n_presence )

lm_2way <- lm(data = full_model_outputs,
              formula = pa_AUC ~ pres_method + bg_method + pres_method*bg_method + n_presence )

lm_2way_re <- lme4::lmer(data = full_model_outputs,
                         formula = pa_AUC ~ pres_method + bg_method + pres_method*bg_method + n_presence + (1|species))

lm_3way_re <- lme4::lmer(data = full_model_outputs,
                         formula = pa_AUC ~ pres_method + bg_method + pres_method*bg_method + pres_method*bg_method*n_presence+n_presence + (1|species))


#Compare model fits
AICcmodavg::aictab(cand.set = list(lm_3way,lm_2way),
                   modnames = c("lm w 3-way int","lm w int"),
                   second.ord = F)

AICcmodavg::aictab(cand.set = list(lm_3way_re,lm_2way_re),
                   modnames = c("re w 3-way int","re w int"),
                   second.ord = F)

summary(lm_2way)
summary(lm_3way)
summary(lm_2way_re)
summary(lm_3way_re)

AIC(lm_2way)
AIC(lm_2way_re)


lsmeans::lsmeans(lm_3way, pairwise~pres_method*bg_method, adjust="tukey")
lsmeans::lsmeans(lm_2way, pairwise~pres_method*bg_method, adjust="tukey")
lsmeans::lsmeans(lm_2way_re, pairwise~pres_method*bg_method, adjust="tukey")

library(nlme)
test<-lme(fixed = pa_AUC ~ pres_method + bg_method + pres_method*bg_method + pres_method*bg_method*n_presence+n_presence,
          random = ~ 1|species,
    data=full_model_outputs, na.action = na.omit)
summary(test)
anova(test)
summary(lm_3way)
summary(lm_2way)
summary(lm_2way_re)


##########################################################

#What determines the relative performance of different bg methods???

group_lookup <- NULL

regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")

for(i in 1:length(regions)){
  
data_i <- disdat::disData(region = regions[i])  
group_lookup <- rbind(group_lookup,
                      unique(data_i$po[c("spid","group")]))
rm(data_i)  
}

#standardize names
group_lookup$group <-
gsub(pattern = "ba",
     replacement = "bats",
     x = group_lookup$group)

group_lookup$group <-
  gsub(pattern = "db",
       replacement = "bird",
       x = group_lookup$group)


group_lookup$group <-
  gsub(pattern = "nb",
       replacement = "bird",
       x = group_lookup$group)

group_lookup$group <-
  gsub(pattern = "ot",
       replacement = "tree",
       x = group_lookup$group)

group_lookup$group <-
  gsub(pattern = "ou",
       replacement = "plant",
       x = group_lookup$group)

group_lookup$group <-
  gsub(pattern = "rt",
       replacement = "tree",
       x = group_lookup$group)

group_lookup$group <-
  gsub(pattern = "ru",
       replacement = "plant",
       x = group_lookup$group)

group_lookup$group <-
  gsub(pattern = "sr",
       replacement = "reptile",
       x = group_lookup$group)

# group_lookup$group <-
#   gsub(pattern = "tree",
#        replacement = "plant",
#        x = group_lookup$group)

length(unique(group_lookup$spid))
unique(group_lookup$group)

full_model_outputs$species

full_model_outputs <-
merge(x = full_model_outputs,
      y= group_lookup,
      by.x = "species",
      by.y="spid")


full_model_outputs %>%
  group_by(pres_method,bg_method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%ungroup%>%
  # mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
  #        pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"),
  #        point_category = fct_relevel(point_category, "low","typical","high")) %>%
  #filter(point_category != "typical")%>%
  ggplot( mapping = aes(x= pa_AUC, fill = group))+
  #ggplot( mapping = aes(x= pa_AUC, fill = pres_method))+
#  geom_histogram(position = "identity", alpha = 0.5)+
  geom_density(alpha=0.5)+
  facet_grid(pres_method ~ bg_method,
             labeller = labeller(
               pres_method = function(x){paste("Presence: \n", x)},
               bg_method = function(x){paste("Background: \n", x)})
  )


###################################################

fold_model_outputs %>%
  mutate(bg_method = fct_relevel(bg_method,"none","rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"))%>%
ggplot(mapping = aes(x=training_AUC,
                     y=testing_AUC))+
  geom_point()+
  facet_grid(pres_method ~ bg_method,
             labeller = labeller(
               pres_method = function(x){paste("Presence: \n", x)},
               bg_method = function(x){paste("Background: \n", x)})
  )+geom_abline(slope = 1,color="blue")+
  xlab("Training AUC")+
  ylab("Testing AUC")

?geom_a

fold_model_outputs$testing_AUC

######################################################

colnames(full_model_output_all)

full_model_output_all %>%
  group_by(method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%ungroup%>%
  # mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
  #        pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"),
  #        point_category = fct_relevel(point_category, "low","typical","high")) %>%
  #filter(point_category != "typical")%>%
  ggplot( mapping = aes(x= pa_AUC,fill=point_category))+
  #ggplot( mapping = aes(x= pa_AUC, fill = pres_method))+
  #  geom_histogram(position = "identity", alpha = 0.5)+
  geom_density(alpha=0.5)+
  facet_wrap("method")
  

ggplot(data = full_model_output_all,
       mapping = aes(x=pa_AUC))+
  geom_density()+
  facet_wrap(facets = "method")



full_model_output_all %>%
  group_by(method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
           ))%>%
  group_by(method,point_category) %>%
  summarise(mean_pa_AUC = mean(pa_AUC),
            median_pa_AUC = median(pa_AUC),
            mean_pa_sensitivity = mean(pa_sensitivity),
            mean_pa_specificity = mean(pa_specificity),
            n_nas = length(which(is.na(pa_AUC)))) -> full_model_summary_by_pa_AUC

write.csv(x = full_model_summary_by_pa_AUC,
          file = "outputs/full_model_summary_by_pa_AUC.csv",
          row.names = F)


colnames(full_model_output_all)





