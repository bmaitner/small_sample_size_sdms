library(disdat)
library(np)
library(AUC)
library(rvinecopulib)
library(pROC)
source("R/get_env_pres.R")
source("R/get_env_bg.R")
source("R/fit_plug_and_play.R")
source("R/pnp_gaussian.R")
source("R/pnp_rangebagging.R")
source("R/pnp_none.R")
source("R/pnp_kde.R")
source("R/project_plug_and_play.R")
source("R/fit_density_ratio.R")
source("R/project_density_ratio.R")
source("R/dr_ulsif.R")
source("R/evaluate_range_map.R")
source("R/get_auc.R")
source("R/pnp_gaussian.R")
source("R/pnp_lobagoc.R")
source("R/pnp_vine.R")
source("R/evaluate_disdat.R")
source("R/stratify_spatial.R")


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
  
  for(i in 1:nrow(models_to_evaluate)){
    
    print("Note: can speed up code considerably by only fitting the background once per location and method,
          since each region uses the same background points")
    
    model_i <- 
    evaluate_disdat(presence_method = models_to_evaluate$presence_method[i],
                    background_method = models_to_evaluate$background_method[i])
    
  
    full_model_outputs <- rbind(full_model_outputs,
                                data.frame(pres_method = models_to_evaluate$presence_method[i],
                                           bg_method = models_to_evaluate$background_method[i],
                                           model_i$full_model_stats))
    
    fold_model_outputs <- rbind(fold_model_outputs,
                                data.frame(pres_method = models_to_evaluate$presence_method[i],
                                           bg_method = models_to_evaluate$background_method[i],
                                           model_i$fold_model_stats))
    
    
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


ggplot(data = full_model_outputs,
       mapping = aes(x= full_AUC, fill = bg_method))+
  geom_histogram()+
  facet_grid(pres_method ~ bg_method)  


##################################
#Plot: pres x bg method PA AUC histograms

colnames(full_model_outputs)


full_model_outputs %>%
  group_by(pres_method,bg_method) %>%
  mutate(min_threshold = stats::quantile(n_presence,.1,na.rm=T),
         max_threshold = stats::quantile(n_presence,.9,na.rm=T))%>%
  mutate(point_category = 
           case_when(n_presence < min_threshold ~ "low",
                     n_presence > max_threshold ~ "high",
                     n_presence >= min_threshold & n_presence <= max_threshold ~ "typical"
                       ))%>%ungroup%>%
  mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"),
         point_category = fct_relevel(point_category, "low","typical","high")) %>%
  #filter(point_category != "typical")%>%
  ggplot( mapping = aes(x= pa_AUC, fill = point_category))+
  #ggplot( mapping = aes(x= pa_AUC, fill = pres_method))+
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
  mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde")) %>%
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


colnames(full_model_outputs)

# full_model_outputs$pres_method <- factor(full_model_outputs$pres_method,levels = c("rangebagging","gaussian","kde"),ordered = T)
# full_model_outputs$bg_method <- factor(full_model_outputs$bg_method,levels = c("rangebagging","gaussian","kde"),ordered = T)
# full_model_outputs$pres_method <- factor(full_model_outputs$pres_method,ordered = F)
# full_model_outputs$bg_method <- factor(full_model_outputs$bg_method,ordered = F)



colnames(full_model_outputs)


summary(lm(data = full_model_outputs,formula = full_AUC ~ pres_method + bg_method + pres_method*bg_method + n_presence))
full_model_outputs$pres_method

colnames(full_model_outputs)

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
  mutate(bg_method = fct_relevel(bg_method,"rangebagging","gaussian","kde"),
         pres_method = fct_relevel(pres_method,"rangebagging","gaussian","kde"))%>%
ggplot(mapping = aes(x=training_AUC,
                     y=testing_AUC))+
  geom_point()+
  facet_grid(pres_method ~ bg_method,
             labeller = labeller(
               pres_method = function(x){paste("Presence: \n", x)},
               bg_method = function(x){paste("Background: \n", x)})
  )+geom_abline(slope = 1,color="blue")

?geom_a

fold_model_outputs$testing_AUC



