# This script is designed to evaluate which of the performance metrics are useful estimates of the PA equivalent metrics

#load packages
library(confintr)


library(tidyverse)
library(ggplot2)
library(ggpubr)
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

tempfile_fold <- "outputs/temp_bakeoff_output_fold.rds"
tempfile_fold_dr <- "outputs/temp_bakeoff_output_fold_dr.rds"

fold_output <- readRDS(tempfile_fold)
fold_output_dr <- readRDS(tempfile_fold_dr)

fold_output_dr %>%
  bind_rows(fold_output) -> fold_output

fold_output %>%
  mutate(model = case_when(!is.na(ratio_method) ~ ratio_method,
                           is.na(ratio_method) ~ paste(pres_method, "/", bg_method))) -> fold_output

rm(fold_output_dr)

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

full_output %>%full_join(mean_fold_output,
                         by = c("species","model") )-> combined_output

#full_output ->combined_output
#########

#How many species have 20 or fewer occurrences?
combined_output %>%
  dplyr::select(species,n_presence)%>%
  na.omit()%>%
  unique()%>%
  dplyr::filter(n_presence <= 20)


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
            #`CI Low Corr. full vs. PA AUC` =  ci_cor(full_AUC,pa_AUC)$interval[1],
            #`CI High Corr. testing vs. PA AUC` = ci_cor(x=full_AUC,y=pa_AUC)$interval[2],
            `mean PA AUC` = mean(pa_AUC,na.rm = TRUE),
            #`CI Low mean PA AUC` = quantile(x = pa_AUC,0.025,na.rm=TRUE),
            #`CI Low mean PA AUC` =  mean(pa_AUC) - (sd(pa_AUC)*1.96),
            #`CI High mean PA AUC` = quantile(x = pa_AUC,0.975,na.rm=TRUE),
            `mean full AUC` = mean(full_AUC,na.rm = TRUE),
            `Corr. full vs. PA AUC` = cor(x=full_AUC,y=pa_AUC,use = "pairwise.complete.obs"),
            cor_full_pa_pred_acc = cor(x = full_prediction_accuracy,y=pa_prediction_accuracy,use = "pairwise.complete.obs"),
            mean_pred_acc_full = mean(full_prediction_accuracy),
            mean_pred_acc_pa = mean(pa_prediction_accuracy),
            cor_full_pa_sensitivity = cor(x = full_sensitivity,y=pa_sensitivity, use = "pairwise.complete.obs"),
            cor_full_pa_specificity = cor(x = full_specificity,y=pa_specificity, use = "pairwise.complete.obs"),
            cor_full_pa_correlation = cor(x = full_correlation,y=pa_correlation, use = "pairwise.complete.obs"),
            cor_full_pa_kappa = cor(x = full_kappa,y=pa_kappa, use = "pairwise.complete.obs")
            ) %>% arrange(`% AUC < 0.50`) -> sss_model_screening

                      
sss_model_screening[2:ncol(sss_model_screening)] <- sss_model_screening[2:ncol(sss_model_screening)] %>% round(digits = 3)

sss_model_screening %>%
  write.csv(file = "tables/sss_model_screening.csv",
            row.names = FALSE)

########################################

# Perhaps a figure? Plot full vs PA for each metric.
# Could facet wrap by metric, color by model? 
# Or models as separate panels, metrics as different figures?

# needed structure
  # model
  # metric
  # pa value
  # full value


colnames(combined_output)


combined_output %>%
  filter(n_presence <= 20)%>%
  group_by(model)%>%
  filter(!all(is.na(pa_AUC)))%>%
  mutate(bad_pa_auc = case_when(pa_AUC < 0.5 ~ 1,
                                pa_AUC >= 0.5 ~ 0))%>%
  pivot_longer(cols = c(full_AUC,pa_AUC),
               names_to = "AUC_type",values_to = "AUC_value") %>%
  pivot_longer(cols = c(full_sensitivity,pa_sensitivity),
               names_to = "sensitivity_type",
               values_to = "sensitivity_value")%>%
  pivot_longer(cols = c(full_specificity, pa_specificity),
               names_to = "specificity_type",
               values_to = "specificity_value")%>%
  pivot_longer(cols = c(full_prediction_accuracy,pa_prediction_accuracy),
             names_to = "pred_acc_type",
             values_to = "pred_acc_value")%>%
  pivot_longer(cols = c(full_correlation,pa_correlation),
               names_to = "correlation_type",
               values_to = "correlation_value")%>%
  pivot_longer(cols = c(full_kappa,pa_kappa),
               names_to = "kappa_type",
               values_to = "kappa_value")->test


  test %>%
  select(c("runtime","entropy","model","AUC_type","AUC_value",
           "sensitivity_type","sensitivity_value",
           "specificity_type","specificity_value",
           "pred_acc_type","pred_acc_value",
           "correlation_type","correlation_value",
           "kappa_type","kappa_value")) %>%
    unique()->test
  

# model
# metric
# pa value
# full value


# AUC

combined_output %>%
  filter(n_presence <= 20) %>%
  group_by(model) %>%
  #filter(!all(is.na(pa_AUC))) %>%
  dplyr::select(model,full_AUC,pa_AUC)%>%
  mutate(metric = "AUC")%>%
  rename(full_value = full_AUC,
         pa_value = pa_AUC)%>%
  
bind_rows(
  combined_output %>%
  filter(n_presence <= 20) %>%
  group_by(model) %>%
  #filter(!all(is.na(pa_AUC))) %>%
  dplyr::select(model,full_sensitivity,pa_sensitivity)%>%
  mutate(metric = "sensitivity")%>%
  rename(full_value = full_sensitivity,
         pa_value = pa_sensitivity))%>%

bind_rows(
  combined_output %>%
    filter(n_presence <= 20) %>%
    group_by(model) %>%
    #filter(!all(is.na(pa_AUC))) %>%
    dplyr::select(model,full_specificity,pa_specificity)%>%
    mutate(metric = "specificity")%>%
    rename(full_value = full_specificity,
           pa_value = pa_specificity))%>%

  bind_rows(
    combined_output %>%
      filter(n_presence <= 20) %>%
      group_by(model) %>%
      #filter(!all(is.na(pa_AUC))) %>%
      dplyr::select(model,full_prediction_accuracy,pa_prediction_accuracy)%>%
      mutate(metric = "prediction_accuracy")%>%
      rename(full_value = full_prediction_accuracy,
             pa_value = pa_prediction_accuracy))%>%
  bind_rows(
    combined_output %>%
      filter(n_presence <= 20) %>%
      group_by(model) %>%
      #filter(!all(is.na(pa_AUC))) %>%
      dplyr::select(model,full_kappa,pa_kappa)%>%
      mutate(metric = "kappa")%>%
      rename(full_value = full_kappa,
             pa_value = pa_kappa)) %>%
  bind_rows(
    combined_output %>%
      filter(n_presence <= 20) %>%
      group_by(model) %>%
      #filter(!all(is.na(pa_AUC))) %>%
      dplyr::select(model,full_correlation,pa_correlation)%>%
      mutate(metric = "correlation")%>%
      rename(full_value = full_correlation,
             pa_value = pa_correlation)) -> metric_corr_data


# Plot of correlations
correlations_training_vs_pa <-
metric_corr_data %>%
  ggplot(mapping = aes(x = full_value,
                       y = pa_value))+
  geom_point(aes(color = model))+
  stat_cor(method = "pearson",)+
  geom_abline(slope = 1,intercept = 0)+
  facet_wrap(~metric)+
  theme_bw()+
  xlab("Value (Training Data)")+
  ylab("Value (P/A Data)")

ggsave(filename = "figures/training_vs_pa_metric_correlations.jpeg",
       plot = correlations_training_vs_pa,
       dpi = 600,width = 10,height = 5,units = "in")


###################################################

# Testing v training correlations

data_for_cor_test<-
combined_output %>%
  filter(n_presence <= 20)%>%
  select(`mean testing_AUC`,pa_AUC)%>%
  na.omit()
  
  
  cor.test(x = data_for_cor_test$`mean testing_AUC`,
           y = data_for_cor_test$pa_AUC,
           method = "pearson")

  
# Checking whether rank orders are correlated  
  
data_for_rank_cor_test <-  
combined_output %>%
  filter(n_presence <= 20) %>%
  select(species,model,`mean testing_AUC`,pa_AUC) %>%
  na.omit()%>%
  group_by(species)%>%
  arrange(species,-pa_AUC)%>%
  mutate(pa_rank = row_number())%>%
  arrange(species,-`mean testing_AUC`)%>%
  mutate(test_rank = row_number())

cor.test(x = data_for_rank_cor_test$test_rank,
         y = data_for_rank_cor_test$pa_rank,
         method = "pearson")

# How often was the model that ranked highest on test data also the best on PA data?

nrow(data_for_rank_cor_test %>%
  filter(test_rank == 1 & pa_rank == 1))/length(unique(data_for_rank_cor_test$species))*100







