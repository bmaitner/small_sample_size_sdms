# The goal of this script is to evaluate the performance of an ensemble of models for small sample size species (which will be compared to the model bakeoff and rarified bakeoff)

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
library(disdat)
source("R/evaluate_disdat_ensemble.R")


# Select models to use in ensemble

  model_vector <- c("maxnet","rulsif","kde/kde")

# Get ensemble performance data

  ensemble_performance <- evaluate_ensemble_disdat(model_vector = model_vector,
                         quantile = 0.05,
                         verbose = TRUE,
                         ncl = 5,
                         temp_file = "outputs/temp_ensemble_eval.RDS")
  
# Get non-ensemble data for comparison
  
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
  
  rm(full_output_dr,tempfile_full,tempfile_full_dr)
  
# Filter non-ensemble stuff to relevant models

  full_output %>%
    mutate(model = gsub(pattern = " ",replacement="",x=model))%>%
    filter(model %in% model_vector) -> relevant_full_output

  
#################################  

# Preliminary visualizations
  
## folds
  ensemble_fold <- ensemble_performance$fold_model_stats
  
# Vote AUCs outperform ensemble AUCs
  ensemble_fold %>%
    ggplot(mapping = aes(x=testing_AUC,
                         y=testing_vote_AUC))+
    geom_point()+
    geom_abline(slope = 1,intercept = 0)

#   
  ensemble_fold %>%
    filter(n_presence <= 20)%>%
    ggplot(mapping = aes(x=testing_AUC,
                         y=testing_vote_AUC))+
    geom_point()+
    geom_abline(slope = 1,intercept = 0)
  
# Pred acc (always higher with ensemble averaging vs vote thresholding)  
  
  ensemble_fold %>%
    ggplot(mapping = aes(x=testing_prediction_accuracy,
                         y=testing_vote_prediction_accuracy))+
    geom_point()+
    geom_abline(slope = 1,intercept = 0)
  
# Sens

  ensemble_fold %>%
    ggplot(mapping = aes(x=testing_sensitivity,
                         y=testing_vote_sensitivity))+
    geom_point()+
    geom_abline(slope = 1,intercept = 0)

# Spec
  
  ensemble_fold %>%
    ggplot(mapping = aes(x=testing_specificity,
                         y=testing_vote_specificity))+
    geom_point()+
    geom_abline(slope = 1,intercept = 0)


###################
  
# Full model
  
  ensemble_full <- ensemble_performance$full_model_stats

  
# Sens
  
# higher sensitivity by thresholding the averaged models
  # at lower densities, seems to be better to take any vote approach

  #species,presence,name(model),value (metric = kde sensitivity)

  
   ensemble_full %>%
     select(species,n_presence, contains("sensitivity")) %>%
     select(species,n_presence, contains("pa_")) %>%
     select(species,n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("sensit"))%>%
     bind_rows(relevant_full_output %>%
                 select(species,n_presence, model, pa_sensitivity) %>%
                 mutate(model = paste(model,"pa_sensitivity")) %>%
                 rename(value = pa_sensitivity,
                        name = model)%>%
                 select(species,n_presence,name,value))%>%
     ggplot(mapping = aes(x = log10(n_presence),
                          y = value,
                          color = name))+
     geom_point()+
     geom_smooth(method = "loess")

   # sensitivity: highest when taking any vote, lowest when taking all votes
   ensemble_full %>%
    filter(n_presence <= 20)%>%
     select(n_presence, contains("sensitivity")) %>%
     select(n_presence, contains("pa_")) %>%
     select(n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("sensit"))%>%
     ggplot(mapping = aes(x = n_presence,
                          y = value,
                          color = name))+
     geom_point()+
     geom_smooth(method = "lm")
   
   
 
# Highest specificity when taking locations that all models agree upon


   ensemble_full %>%
     select(n_presence, contains("specificity")) %>%
     select(n_presence, contains("pa_")) %>%
     select(n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("spec"))%>%
     ggplot(mapping = aes(x = log10(n_presence),
                          y = value,
                          color = name))+
     geom_point()+
     geom_smooth(method = "loess")

  # for low sample sizes, vote ensembles are all above .9, where averaging can be quite poor (0.3)  
   
   ensemble_full %>%
     filter(n_presence <= 20)%>%
     select(n_presence, contains("specificity")) %>%
     select(n_presence, contains("pa_")) %>%
     select(n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("spec"))%>%
     ggplot(mapping = aes(x = n_presence,
                          y = value,
                          color = name))+
     geom_point()+
     geom_smooth(method = "lm")
   
   
# Method pred accuracy vs sample size - best when taking where all models agree
   
   ensemble_full %>%
     select(n_presence, contains("prediction_acc")) %>%
     select(n_presence, contains("pa_")) %>%
     pivot_longer(cols = contains("pred"))%>%
     ggplot(mapping = aes(x = log10(n_presence),
                          y = value,
                          color = name))+
     geom_point()+
     geom_smooth(method = "loess")
   
   ensemble_full %>%
     filter(n_presence <= 20) %>%
     select(n_presence, contains("prediction_acc")) %>%
     select(n_presence, contains("pa_")) %>%
     pivot_longer(cols = contains("pred"))%>%
     ggplot(mapping = aes(x = n_presence,
                          y = value,
                          color = name))+
     geom_point()+
     geom_smooth(method = "loess")

   
# Plot: bar plot of prediction accuracy, sensitivity, specificity with different votes and component methods

   
   ensemble_full %>%
     select(species,n_presence, contains("sensitivity")|contains("specificity")|contains("prediction")) %>%
     select(species,n_presence, contains("pa_")) %>%
     select(species,n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("sensit")|contains("specif")|contains("prediction"))%>%
     mutate(model = case_when(grepl(pattern = "any_vote",x=name) ~ "Ensemble_any_support",
                              grepl(pattern = "all_vote",x=name) ~ "Ensemble_unanimous_support",
                              .default = "Ensemble_mean"))%>%
     mutate(metric = case_when(grepl(pattern = "sensitivity",x=name) ~ "sensitivity",
                               grepl(pattern = "specificity",x=name) ~ "specificity",
                               grepl(pattern = "prediction",x=name) ~ "prediction_accuracy"))%>%
     bind_rows(   relevant_full_output %>%
                    select(species,n_presence, model, pa_sensitivity,pa_specificity,pa_prediction_accuracy) %>%
                    pivot_longer(contains("pa_"),names_to = "metric")%>%
                    mutate(metric = gsub(pattern = "pa_",replacement = "",x=metric))) %>%
     # mutate(model = factor(x=model,levels = c("ensemble_all_votes","ensemble_mean","ensemble_any_votes",
     #                                          "kde/kde","rulsif","maxnet")))%>%
     mutate(metric = case_when(grepl(pattern = "sensitivity",x=metric) ~ "Sensitivity",
                               grepl(pattern = "specificity",x=metric) ~ "Specificity",
                               grepl(pattern = "prediction",x=metric) ~ "Prediction accuracy"))%>%
     mutate(model = gsub(pattern = "_",replacement =" ",x=model))%>%
     mutate(metric = gsub(pattern = "_",replacement =" ",x=metric))%>%
     mutate(model = case_when(grepl(pattern = "rulsif",x=model) ~ "ruLSIF",
                               grepl(pattern = "maxnet",x=model) ~ "Maxnet",
                               grepl(pattern = "kde/kde",x=model) ~ "KDE/KDE",
                              .default = model))%>%
     mutate(model = factor(x=model,levels = c("Ensemble unanimous support","KDE/KDE",
                                              "ruLSIF","Ensemble mean",
                                              "Maxnet","Ensemble any support")))%>%
     mutate(metric = factor(x=metric,
                            levels = c("Sensitivity","Specificity","Prediction accuracy")))%>%
    filter(n_presence <= 20)%>%
   ggplot(mapping = aes(x = model,
                          y = value))+
     geom_hline(data = . %>%
                  group_by(metric,model)%>%
                  summarize(med_value = median(value))%>%
                  slice_max(order_by = med_value,n = 1),
                mapping = aes(yintercept = med_value),lty=2)+
     geom_boxplot()+
     facet_wrap(~metric,ncol = 1)+
     theme_bw()+
     xlab(NULL)+
     ylab(NULL)+
     theme(strip.text = element_text(size = 10 ),
           axis.text.x = element_text(size = 10 )) -> ensemble_figure
   
ensemble_figure     

ggsave(filename = "figures/ensemble_metrics.jpg",
       plot = ensemble_figure,width = 10,height = 5,units = "in",dpi = 600)
   
      
# Need a table equivalent to table 3 for ensembles
   
   
   ensemble_full %>%
     select(species,n_presence, contains("sensitivity")|
              contains("specificity")|
              contains("prediction")|
              contains("correlation")|
              contains("kappa")) %>%
     select(species,n_presence, contains("pa_")) %>%
     select(species,n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("pa_"))%>%
     mutate(model = case_when(grepl(pattern = "any_vote",x=name) ~ "ensemble_any_votes",
                              grepl(pattern = "all_vote",x=name) ~ "ensemble_all_votes",
                              .default = "ensemble_mean"))%>%
     mutate(metric = case_when(grepl(pattern = "sensitivity",x=name) ~ "sensitivity",
                               grepl(pattern = "specificity",x=name) ~ "specificity",
                               grepl(pattern = "prediction",x=name) ~ "prediction_accuracy",
                               grepl(pattern = "correlation",x=name) ~ "correlation",
                               grepl(pattern = "kappa",x=name) ~ "kappa"))%>%
     select(species,n_presence,value,model,metric)%>%
     bind_rows(   relevant_full_output %>%
                    select(species,n_presence, model,
                           pa_sensitivity,
                           pa_specificity,
                           pa_prediction_accuracy,
                           pa_correlation,
                           pa_kappa) %>%
                    pivot_longer(contains("pa_"),names_to = "metric")%>%
                    mutate(metric = gsub(pattern = "pa_",replacement = "",x=metric))) %>%
     # mutate(model = factor(x=model,levels = c("ensemble_all_votes","ensemble_mean","ensemble_any_votes",
     #                                          "kde/kde","rulsif","maxnet")))%>%
     mutate(model = gsub(pattern = "_",replacement =" ",x=model))%>%
     mutate(metric = gsub(pattern = "_",replacement =" ",x=metric))%>%
     mutate(model = factor(x=model,levels = c("ensemble all votes","kde/kde",
                                              "ensemble mean","rulsif",
                                              "maxnet","ensemble any votes")))%>%
     #mutate(metric = factor(x=metric,levels = c("sensitivity","specificity","prediction accuracy")))%>%
     filter(n_presence <= 20)%>%
     filter(metric != "correlation")%>%
     group_by(model,metric) %>%
     summarise(mean_value = mean(value))%>%
     pivot_wider(names_from = "metric", values_from = "mean_value")%>%
     arrange(-specificity)%>%
     dplyr::select(model,`prediction accuracy`,specificity,sensitivity,kappa)%>%
     mutate(across(1:4, function(x){round(x,digits = 3)}))-> ensemble_table
   
   ensemble_table
   

   write.csv(x = ensemble_table,
             file = "tables/ensemble_table.csv",
             row.names = FALSE)   

   
#############################################

# Ensemble spanning less sensitivity-specificity  
   
   # lets keep maxnet, but shrink the other side
   # my guess is that, since these are all on the sensitivity side, 
    # the ensemble will perform less well on specificity, but comparably on sensitivity
   
   less_variation_ensemble_performance <- evaluate_ensemble_disdat(model_vector = c("maxnet","gaussian/gaussian","gaussian/kde"),
                                                                   quantile = 0.05,
                                                                   verbose = TRUE,
                                                                   ncl = 5,
                                                                   temp_file = "outputs/temp_low_variation_ensemble_eval.RDS")
   
   less_variation_ensemble_performance$full_model_stats %>%
     select(species,n_presence, contains("sensitivity")|contains("specificity")|contains("prediction")) %>%
     select(species,n_presence, contains("pa_")) %>%
     select(species,n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("sensit")|contains("specif")|contains("prediction"))%>%
     mutate(model = case_when(grepl(pattern = "any_vote",x=name) ~ "ensemble_any_votes",
                              grepl(pattern = "all_vote",x=name) ~ "ensemble_all_votes",
                              .default = "ensemble_mean"))%>%
     mutate(metric = case_when(grepl(pattern = "sensitivity",x=name) ~ "sensitivity",
                               grepl(pattern = "specificity",x=name) ~ "specificity",
                               grepl(pattern = "prediction",x=name) ~ "prediction_accuracy"))%>%
     bind_rows(     full_output %>%
                      mutate(model = gsub(pattern = " ",replacement="",x=model))%>%
                      filter(model %in% c("maxnet","gaussian/gaussian","gaussian/kde"))%>%
                    select(species,n_presence, model, pa_sensitivity,pa_specificity,pa_prediction_accuracy) %>%
                    pivot_longer(contains("pa_"),names_to = "metric")%>%
                    mutate(metric = gsub(pattern = "pa_",replacement = "",x=metric))) %>%
     # mutate(model = factor(x=model,levels = c("ensemble_all_votes","ensemble_mean","ensemble_any_votes",
     #                                          "kde/kde","rulsif","maxnet")))%>%
     mutate(model = gsub(pattern = "_",replacement =" ",x=model))%>%
     mutate(metric = gsub(pattern = "_",replacement =" ",x=metric))%>%
     # mutate(model = factor(x=model,levels = c("ensemble all votes","kde/kde",
     #                                          "ensemble mean","rulsif",
     #                                          "maxnet","ensemble any votes")))%>%
     mutate(metric = factor(x=metric,levels = c("sensitivity","specificity","prediction accuracy")))%>%
     filter(n_presence <= 20)%>%
     ggplot(mapping = aes(x = model,
                          y = value))+
     geom_hline(data = . %>%
                  group_by(metric,model)%>%
                  summarize(med_value = median(value))%>%
                  slice_max(order_by = med_value,n = 1),
                mapping = aes(yintercept = med_value),lty=2)+
     geom_boxplot()+
     facet_wrap(~metric,ncol = 1)+
     theme_bw()+
     xlab(NULL)+
     ylab(NULL)
   
   
   less_variation_ensemble_performance$full_model_stats %>%
     select(species,n_presence, contains("sensitivity")|
              contains("specificity")|
              contains("prediction")|
              contains("correlation")|
              contains("kappa")) %>%
     select(species,n_presence, contains("pa_")) %>%
     select(species,n_presence, !contains("pAUC_")) %>%
     pivot_longer(cols = contains("pa_"))%>%
     mutate(model = case_when(grepl(pattern = "any_vote",x=name) ~ "ensemble_any_votes",
                              grepl(pattern = "all_vote",x=name) ~ "ensemble_all_votes",
                              .default = "ensemble_mean"))%>%
     mutate(metric = case_when(grepl(pattern = "sensitivity",x=name) ~ "sensitivity",
                               grepl(pattern = "specificity",x=name) ~ "specificity",
                               grepl(pattern = "prediction",x=name) ~ "prediction_accuracy",
                               grepl(pattern = "correlation",x=name) ~ "correlation",
                               grepl(pattern = "kappa",x=name) ~ "kappa"))%>%
     select(species,n_presence,value,model,metric)%>%
     # bind_rows(   relevant_full_output %>%
     #                select(species,n_presence, model,
     #                       pa_sensitivity,
     #                       pa_specificity,
     #                       pa_prediction_accuracy,
     #                       pa_correlation,
     #                       pa_kappa) %>%
     #                pivot_longer(contains("pa_"),names_to = "metric")%>%
     #                mutate(metric = gsub(pattern = "pa_",replacement = "",x=metric))) %>%
     # mutate(model = factor(x=model,levels = c("ensemble_all_votes","ensemble_mean","ensemble_any_votes",
     #                                          "kde/kde","rulsif","maxnet")))%>%
     mutate(model = gsub(pattern = "_",replacement =" ",x=model))%>%
     mutate(metric = gsub(pattern = "_",replacement =" ",x=metric))%>%
     mutate(model = factor(x=model,levels = c("ensemble all votes","kde/kde",
                                              "ensemble mean","rulsif",
                                              "maxnet","ensemble any votes")))%>%
     #mutate(metric = factor(x=metric,levels = c("sensitivity","specificity","prediction accuracy")))%>%
     filter(n_presence <= 20)%>%
     filter(metric != "correlation")%>%
     group_by(model,metric) %>%
     summarise(mean_value = mean(value))%>%
     pivot_wider(names_from = "metric", values_from = "mean_value")%>%
     arrange(-specificity)%>%
     dplyr::select(model,`prediction accuracy`,specificity,sensitivity,kappa)%>%
     mutate(across(1:4, function(x){round(x,digits = 3)}))-> ensemble_table_less_variation
   
##############################################
   
# Ensemble spanning same sens-specificity, but with more models
   
   # selected two additional models that were intermediate between each of the existing models for a total of 7
   
   # my guess is that the all/any votes should be the same, but with longer model runtimes
   
   # Get ensemble performance data
   
   lots_of_models_ensemble_performance <- evaluate_ensemble_disdat(model_vector = c("maxnet",
                                                                     "gaussian/gaussian",
                                                                     "rangebagging/none",
                                                                     "rulsif",
                                                                     "lobagoc/none",
                                                                     "ulsif",
                                                                     "kde/kde"),
                                                    quantile = 0.05,
                                                    verbose = TRUE,
                                                    ncl = 5,
                                                    temp_file = "outputs/temp_lots_of_models_ensemble_eval.RDS")
   
#######################################################
   library(grid)
   
    poor_models <- read.csv("tables/small_sample_size_comparison_to_maxnet.csv") %>%
    filter(is.na(pval) | pval <= 0.05)%>%
     pull(model)
   
      # combine ensemble outputs
   
    # should be sensitivity with any votes, specificity with all votes
   
   
     ensemble_performance$full_model_stats %>%
     select(species,n_presence,
            pa_any_vote_sensitivity,
            pa_any_vote_specificity,
            pa_all_vote_specificity,
            pa_all_vote_sensitivity) %>%
     mutate(model = "ensemble_wide_few")%>%
     bind_rows( less_variation_ensemble_performance$full_model_stats %>%
                  select(species,n_presence,
                         pa_any_vote_sensitivity,
                         pa_any_vote_specificity,
                         pa_all_vote_specificity,
                         pa_all_vote_sensitivity) %>%
                  mutate(model = "ensemble_narrow_few"))%>%
     bind_rows( lots_of_models_ensemble_performance$full_model_stats %>%
                  select(species,n_presence,
                         pa_any_vote_sensitivity,
                         pa_any_vote_specificity,
                         pa_all_vote_specificity,
                         pa_all_vote_sensitivity) %>%
                  mutate(model = "ensemble_wide_many"))%>%
     mutate(All = pa_all_vote_sensitivity-pa_all_vote_specificity,
            Any = pa_any_vote_sensitivity-pa_any_vote_specificity,
            All_sum = pa_all_vote_sensitivity + pa_all_vote_specificity,
            Any_sum = pa_any_vote_sensitivity + pa_any_vote_specificity)%>%
       pivot_longer(cols = c(All,Any),
                    names_to = "Votes",
                    values_to = "sens_spec_ratio") %>%
       pivot_longer(cols = c(All_sum,Any_sum),
                    names_to = "Votes_sum",
                    values_to = "sens_spec_sum") %>%
       mutate(Votes_sum = gsub(pattern = "_sum",
                               replacement = "",
                               x=Votes_sum)) %>%
       filter(Votes == Votes_sum) %>%
       # join with the non-ensemble data
       
       bind_rows( full_output %>%
                          select(species,n_presence,
                                 pa_sensitivity,
                                 pa_specificity,
                                 model)%>%
                          mutate(Votes="NA",
                                 sens_spec_ratio = pa_sensitivity - pa_specificity,
                                 sens_spec_sum = pa_sensitivity + pa_specificity)
                        ) %>%
       #filter to rare species
       filter(n_presence <=20) %>%
       #filter out anything that did worse than maxent for rare species
       filter(!model %in% poor_models)%>%
       filter(model != "CVmaxnet")%>%
       group_by(model,Votes) %>%
       summarise(mean_sens_spec = median(sens_spec_ratio,na.rm=TRUE),
                 sens_spec_low = quantile(sens_spec_ratio,0.25,na.rm=TRUE),
                 sens_spec_high = quantile(sens_spec_ratio,0.75,na.rm=TRUE),
                 med_sum_sens_spec = median(sens_spec_sum,na.rm=TRUE) ) -> temp_data
     

     temp_data%>%
     ggplot(mapping = aes(y=model,x=mean_sens_spec,color=Votes))+
     geom_point(mapping = aes(size=med_sum_sens_spec))+
      geom_errorbarh(mapping = aes(y=model,
                                   xmin=sens_spec_low,
                                   xmax=sens_spec_high,
                                   color=Votes))+
       xlim(c(-1,1))+
       scale_y_discrete(limits=temp_data%>%
                          group_by(model)%>%
                          summarize(order = case_when(all(Votes == 'NA')~ mean(mean_sens_spec),
                                                      all(Votes != "NA") ~ 1+var(mean_sens_spec)))%>%
                          arrange(-order)%>%
                          pull(model))+
       xlab("Sensitivity - Specificity")+
       labs(size = "Sensitivity +\nSpecificity") -> temp
     temp 
     
     #text_high <- textGrob("High Specificity", gp=gpar(fontsize=13, fontface="bold"))
     
     text_high_spec <- textGrob("High Specificity\nPresences Correct\nAssumes good sampling")
     text_high_sens <- textGrob("High Sensitivity\nAbsences Correct\nAssumes poor sampling")
     

     
     
     temp <- temp+
       annotation_custom(text_high_spec,xmin=-1,xmax=-1,ymin=-8,ymax=-6)+
       annotation_custom(text_high_sens,xmin=1,xmax=1,ymin=-8,ymax=-6)+
       theme(plot.margin = unit(c(1,1,4,1), "lines")) + #top,right,bottom,left
       coord_cartesian(ylim=c(25,0), clip="off")
     
     
    ggsave(filename = "figures/ensemble_sens_spec.jpg",
           plot = temp,
           width = 10,
           height = 4.5,
           units = "in",
           dpi = 600)
    
    ggsave(filename = "figures/ensemble_sens_spec.svg",
           plot = temp,
           width = 10,
           height = 4.5,
           units = "in",
           dpi = 600)   
    

     

     
     