#load packages
library(confintr)
# Load libraries

library(tidyverse)
library(ggplot2)
library(tidyverse)

# Load and format full data

tempfile_full <- "outputs/temp_bakeoff_output_full.rds"
tempfile_full_dr <- "outputs/temp_bakeoff_output_full_dr.rds"


full_output <- readRDS(tempfile_full)
full_output_dr <- readRDS(tempfile_full_dr)

full_output %>%
  mutate(pa_AUC = as.numeric(pa_AUC),
         full_AUC = as.numeric(full_AUC)) -> full_output

full_output_dr %>%
  mutate(pa_AUC = as.numeric(pa_AUC),
         full_AUC = as.numeric(full_AUC)) -> full_output_dr

full_output_dr %>%
  bind_rows(full_output)->full_output

full_output %>%
  mutate(model = case_when(!is.na(ratio_method) ~ ratio_method,
                           is.na(ratio_method) ~ paste(pres_method, "/", bg_method))) -> full_output

rm(full_output_dr)

# combining fold data with full data

fold_output %>%
  group_by(model,species)%>%
  summarise('mean training AUC' = na.omit(training_AUC) %>% mean(),
            'mean training_pAUC_specificity' = na.omit(training_pAUC_specificity) %>% mean(),
            'mean training_pAUC_sensitivity' = na.omit(training_pAUC_sensitivity) %>% mean(),
            'mean testing_pAUC_specificity' = na.omit(testing_pAUC_specificity) %>% mean(),
            'mean testing_pAUC_sensitivity' = na.omit(testing_pAUC_sensitivity) %>% mean(),
            'mean testing_DOR' = na.omit(testing_DOR) %>% mean(),
            'mean testing_prediction_accuracy' = na.omit(testing_prediction_accuracy) %>% mean(),
            'mean testing_AUC' = na.omit(testing_AUC) %>% mean(),
            'mean testing_sensitivity' = na.omit(testing_sensitivity) %>% mean(),
            'mean testing_specificity' = na.omit(testing_specificity) %>% mean(),
            'mean testing_correlation' = na.omit(testing_correlation) %>% mean(),
            'mean testing_kappa' = na.omit(testing_kappa) %>% mean(),
            'mean entropy' = na.omit(entropy) %>% mean())-> mean_fold_output

full_output %>%
  full_join(mean_fold_output,by = c("species","model"))->combined_output

#########



#################

combined_output %>%
  filter(n_presence <= 20)%>%
  group_by(model)%>%
  filter(all(!is.na(pa_AUC)))%>%
  mutate(bad_pa_auc = case_when(pa_AUC < 0.5 ~ 1,
                                pa_AUC >= 0.5 ~ 0)) %>%
  summarise(`% AUC < 0.50` = sum(bad_pa_auc)/n()*100,
            `Corr. testing vs. PA AUC` = cor(x=`mean testing_AUC`,y=pa_AUC,use = "pairwise.complete.obs"),
            `CI Low Corr. testing vs. PA AUC` =  ci_cor(`mean testing_AUC`,pa_AUC)$interval[1],
            `CI High Corr. testing vs. PA AUC` = ci_cor(x=`mean testing_AUC`,y=pa_AUC)$interval[2],
            `mean PA AUC` = mean(pa_AUC),
            #`CI Low mean PA AUC` = quantile(x = pa_AUC,0.025,na.rm=TRUE),
            `CI Low mean PA AUC` =  mean(pa_AUC) - (sd(pa_AUC)*1.96),
            #`CI High mean PA AUC` = quantile(x = pa_AUC,0.975,na.rm=TRUE),
            `mean testing AUC` = mean(`mean testing_AUC`,na.rm = TRUE)
            ) %>% arrange(`% AUC < 0.50`) -> sss_model_screening

library(plotrix)

sss_model_screening[2:ncol(sss_model_screening)] <- sss_model_screening[2:ncol(sss_model_screening)] %>% round(digits = 3)

sss_model_screening %>%
  write.csv(file = "tables/sss_model_screening.csv",
            row.names = FALSE)

combined_output %>%
  filter(n_presence <= 20)%>%
  group_by(model)

  out <- cor(x=combined_output$full_AUC,y=combined_output$pa_AUC,use="pairwise.complete.obs")
  out <- ci_cor(x=combined_output$full_AUC,y=combined_output$pa_AUC,use="pairwise.complete.obs")$interval[1]
  
out$interval[1] 

  
combined_output %>%
  filter(n_presence <= 20)%>%
  filter(model == "gaussian / lobagoc")%>%
  select(pa_AUC)

quantile(x = combined_output$pa_AUC,0.025,na.rm=TRUE)

#############################################

combined_output %>%
  filter(n_presence <= 20)%>%
  group_by(model)%>%
  filter(!all(is.na(pa_AUC)))%>%
  mutate(bad_pa_auc = case_when(pa_AUC < 0.5 ~ 1,
                                pa_AUC >= 0.5 ~ 0)) %>%
  summarise(`% AUC < 0.50` = sum(bad_pa_auc,na.rm = TRUE)/n()*100,
            `Corr. full vs. PA AUC` = cor(x=full_AUC,y=pa_AUC,use = "pairwise.complete.obs"),
            #`CI Low Corr. full vs. PA AUC` =  ci_cor(full_AUC,pa_AUC)$interval[1],
            #`CI High Corr. testing vs. PA AUC` = ci_cor(x=full_AUC,y=pa_AUC)$interval[2],
            `mean PA AUC` = mean(pa_AUC,na.rm = TRUE),
            #`CI Low mean PA AUC` = quantile(x = pa_AUC,0.025,na.rm=TRUE),
            #`CI Low mean PA AUC` =  mean(pa_AUC) - (sd(pa_AUC)*1.96),
            #`CI High mean PA AUC` = quantile(x = pa_AUC,0.975,na.rm=TRUE),
            `mean full AUC` = mean(full_AUC,na.rm = TRUE),
            #mean_pred_acc_full = mean(full_prediction_accuracy),
            #mean_pred_acc_pa = mean(pa_prediction_accuracy),
            cor_full_pa_pred_acc = cor(x = full_prediction_accuracy,y=pa_prediction_accuracy,use = "pairwise.complete.obs"),
            cor_full_pa_sensitivity = cor(x = full_sensitivity,y=pa_sensitivity, use = "pairwise.complete.obs"),
            cor_full_pa_specificity = cor(x = full_specificity,y=pa_specificity, use = "pairwise.complete.obs"),
            cor_full_pa_correlation = cor(x = full_correlation,y=pa_correlation, use = "pairwise.complete.obs"),
            cor_full_pa_kappa = cor(x = full_kappa,y=pa_kappa, use = "pairwise.complete.obs")
            ) %>% arrange(`% AUC < 0.50`) -> sss_model_screening

                      

##################

# entropy
full_output %>%
  filter(bg_method == "none" | is.na(bg_method))%>%
  filter(n_presence > 20)%>%
  ggplot(mapping = aes(x=log10(n_presence),
                       y=entropy,
                       color=model))+
  geom_point()+
  geom_smooth(method = "lm")

# ent vs spec

  full_output %>%
    filter(bg_method == "none" | is.na(bg_method))%>%
    filter(n_presence > 20)%>%
    ggplot(mapping = aes(x=full_specificity,
                         y=entropy,
                         color=model))+
    geom_point()+
    geom_smooth(method = "lm")

  full_output %>%
    filter(bg_method == "none" | is.na(bg_method))%>%
    #filter(n_presence > 20)%>%
    ggplot(mapping = aes(x=pa_specificity,
                         y=entropy,
                         color=model))+
    geom_point()+
    geom_smooth(method = "lm")

  full_output %>%
    filter(bg_method == "none" | is.na(bg_method))%>%
    filter(n_presence > 20)%>%
    ggplot(mapping = aes(x=pa_sensitivity,
                         y=entropy,
                         color=model))+
    geom_point()+
    geom_smooth(method = "lm")
  
  
  full_output %>%
    filter(bg_method == "none" | is.na(bg_method))%>%
    filter(n_presence > 20)%>%
    ggplot(mapping = aes(x=pa_sensitivity,
                         y=pa_specificity,
                         color=model))+
    geom_point()+
    geom_smooth(method = "lm")

  
  ci_min =0.25
  ci_max = 0.75
  
  full_output %>%
    filter(bg_method == "none" | is.na(bg_method))%>%
    #filter(n_presence > 20)%>%
    group_by(model)%>%
    summarise(mean_pa_sens = mean(pa_sensitivity,na.rm = TRUE),
              ci_high_pa_sens = quantile(pa_sensitivity,probs = ci_max,na.rm = TRUE),
              ci_low_pa_sens = quantile(pa_sensitivity,probs = ci_min,na.rm = TRUE),
              mean_pa_spec = mean(pa_specificity,na.rm = TRUE),
              ci_high_pa_spec = quantile(pa_specificity,probs = ci_max,na.rm = TRUE),
              ci_low_pa_spec = quantile(pa_specificity,probs = ci_min,na.rm = TRUE))%>%
    ggplot(mapping = aes(x=mean_pa_sens,y=mean_pa_spec,color=model))+
    geom_point()+
    geom_errorbar(mapping = aes(ymin=ci_low_pa_spec,
                                ymax = ci_high_pa_spec))+
    geom_errorbarh(mapping = aes(xmin=ci_low_pa_sens,
                                xmax = ci_high_pa_sens))+
    xlab("Sensitivty (P/A)")+
    ylab("Specificity (P/A)")+
    theme_bw()
    
#################
  
  #how many times does each model perform best?
  
    full_output %>%
      select(species,n_presence,model,pa_AUC)%>%
      group_by(species)%>%
      slice(which.max(pa_AUC))%>%
      ungroup()%>%
      group_by(model)%>%
      summarise(n = n())%>%
      arrange(-n)%>%
    mutate(prop_pct = (n/sum(n)) %>% round(2)*100)
  
  full_output %>%
    filter(n_presence <= 20) %>%
    select(species,n_presence,model,pa_AUC)%>%
    group_by(species)%>%
    slice(which.max(pa_AUC))%>%
    ungroup()%>%
    group_by(model)%>%
    summarise(n = n())%>%
    arrange(-n)%>%
    mutate(prop_pct = (n/sum(n)) %>% round(2)*100)

    

    