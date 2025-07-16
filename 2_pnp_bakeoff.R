# 2 pnp bakeoff analyses

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
    bind_rows(full_output) -> full_output
  
  full_output %>%
    mutate(model = case_when(!is.na(ratio_method) ~ ratio_method,
                             is.na(ratio_method) ~ paste(pres_method, "/", bg_method))) -> full_output

  rm(full_output_dr)

# # Load and format cv data

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
  
  
  
#############
  
  # Generate result summaries
  
  
  # Table 1: performance ranked with Presence/Abscence data
  
    Table1 <-
    full_output %>%
      filter(model != "CVmaxnet" ) %>%
      mutate(pa_TSS = pa_sensitivity + pa_specificity - 1) %>%
      group_by(model) %>%
      summarise('median PA AUC' = na.omit(pa_AUC) %>% median(),
                'mean PA AUC' = na.omit(pa_AUC) %>% mean(),
                'mean PA sensitivity' = na.omit(pa_sensitivity) %>% mean(),
                'mean PA specificity' = na.omit(pa_specificity) %>% mean(),
                'mean PA pAUC (sensitivity 0.8 - 1)' = na.omit(pa_pAUC_sensitivity) %>% mean(),
                'mean PA pAUC (specificity 0.8 -1)' = na.omit(pa_pAUC_specificity) %>% mean(),
                'mean PA prediction accuracy' = na.omit(pa_prediction_accuracy) %>% mean(),
                'mean PA correlation' = na.omit(pa_correlation) %>% mean(),
                'mean PA kappa' = na.omit(pa_kappa) %>% mean(),
                'mean PA TSS' = na.omit(pa_TSS) %>% mean()) %>%
    arrange(-`mean PA AUC`)

  Table1[2:ncol(Table1)] <- Table1[2:ncol(Table1)] %>% round(digits = 3)
  
  Table1 %>%
  write.csv(file = "tables/Table1.csv",
            row.names = FALSE)
  
  
##################################  

  # Table 2: as table 1, but limited to small sample sizes only
  Table2 <-
  full_output %>%
    filter(model != "CVmaxnet" ) %>%
    filter(n_presence <= 20) %>%
    mutate(pa_TSS = pa_sensitivity + pa_specificity - 1) %>%
    group_by(model) %>%
    summarise('median PA AUC' = na.omit(pa_AUC) %>% median(),
              'mean PA AUC' = na.omit(pa_AUC) %>% mean(),
              'mean PA sensitivity' = na.omit(pa_sensitivity) %>% mean(),
              'mean PA specificity' = na.omit(pa_specificity) %>% mean(),
              'mean PA pAUC (sensitivity 0.8 - 1)' = na.omit(pa_pAUC_sensitivity) %>% mean(),
              'mean PA pAUC (specificity 0.8 -1)' = na.omit(pa_pAUC_specificity) %>% mean(),
              'mean PA prediction accuracy' = na.omit(pa_prediction_accuracy) %>% mean(),
              'mean PA correlation' = na.omit(pa_correlation) %>% mean(),
              'mean PA kappa' = na.omit(pa_kappa) %>% mean(),
              'mean PA TSS' = na.omit(pa_TSS) %>% mean()) %>%
    arrange(-`mean PA AUC`)
  
  Table2[2:ncol(Table2)] <- Table2[2:ncol(Table2)] %>% round(digits = 3)
  
  Table2 %>%
    write.csv(file = "tables/Table2.csv",
              row.names = FALSE)
  
######################################################################
  
      
        
##############
  
      
      # Mann-Whitney test for comparison of aucs with maxnet (doesn't require a gaussian dist)
      
      for(i in 1:length(unique(full_output$model))){
        
        model_i <- unique(full_output$model)[i]
        
        
        data_x <- full_output %>%
          filter(model == model_i) %>%
          pull(pa_AUC)
        
        data_y <- full_output %>%
          filter(model=="maxnet") %>%
          pull(pa_AUC)
        
        ti <-tryCatch(wilcox.test(x = data_x,
                                  y = data_y),
                      error=function(e){e})
        
        
        if(inherits(ti,"error")){
          ti$p.value <- NA
          ti$statistic <- NA
        }
        
        out_i <- data.frame(model = model_i,
                            comparison = "maxnet",
                            pval = ti$p.value,
                            W = ti$statistic,
                            comparison_w_maxnet  = mean(data_x,na.rm = TRUE) - mean(data_y,na.rm = TRUE))
        
        
        
        if(i==1){ttest_out <- out_i}else(ttest_out <- bind_rows(ttest_out,out_i))
        
        rm(out_i)
        
      }
      
      ttest_out %>%
        mutate(signif = pval <= 0.05) -> models_v_maxent
      
      models_v_maxent %>%
        filter(!signif)%>%
        pull(model) -> all_sample_comparable_to_maxnet
      
      # small sample species t-test vs maxent
      
      for(i in 1:length(unique(full_output$model))){
        
        model_i <- unique(full_output$model)[i]
        
        data_x <- full_output %>%
          filter(n_presence <= 20) %>%
          filter(model==model_i)%>%
          pull(pa_AUC)
        
        data_y <- full_output %>%
          filter(n_presence <= 20) %>%
          filter(model=="maxnet")%>%
          pull(pa_AUC)
        
        ti <-tryCatch(wilcox.test(x = data_x,
                                  y = data_y),
                      error = function(e){e})
        
        if(inherits(ti,"error")){
          ti$p.value <- NA
          ti$statistic <- NA
          
        }
        
        out_i <- data.frame(model = model_i,
                            comparison = "maxnet",
                            pval = ti$p.value,
                            W = ti$statistic,
                            comparison_w_maxnet  = mean(data_x,na.rm = TRUE) - mean(data_y,na.rm = TRUE))
        
        if(i==1){ttest_small_out <- out_i}else(ttest_small_out <- bind_rows(ttest_small_out,out_i))
        
        
        
        rm(out_i)
        
      }
      
      ttest_small_out %>%
        mutate(signif = pval <= 0.05) -> small_samples_models_v_maxent
      
      small_samples_models_v_maxent %>%
        filter(!signif) %>%
        pull(model) -> small_sample_comparable_to_maxnet
      
      small_samples_models_v_maxent %>%
        filter(model != "CVmaxnet" ) %>%
        select(model,W,pval) %>%
        arrange(-pval)%>%
        mutate(pval = round(pval,3)) %>%
        write.csv(file = "tables/small_sample_size_comparison_to_maxnet.csv",
                  row.names = FALSE)
    
      models_v_maxent %>%
        filter(model != "CVmaxnet" ) %>%
        select(model,W,pval) %>%
        arrange(-pval)%>%
        mutate(pval = round(pval,3)) %>%
        write.csv(file = "tables/all_sample_size_comparison_to_maxnet.csv",
                  row.names = FALSE)
      
      # What proportion of algorithms were not significantly worse than maxent?
      
      small_samples_models_v_maxent %>%
        filter(model != "maxnet")%>%
        mutate(signif = as.numeric(signif))%>%
        summarise(Pct_ns_from_maxent = sum(!signif)/n()*100)
        
##################
      
  # which models are significantly different from the best performing in each class?
      
  source("R/sig_different_by_column.R")    
      
      full_output %>%
        select(-n_presence,-n_background,-bg_method,-pres_method,-ratio_method,
               -n_pa_absence,-n_pa_presence,-species)%>%
        select(model,contains("pa_"))%>%
        sig_different_by_column(filter_by = "pa_AUC") -> all_sample_comparable_to_best
      
  write.csv(x = all_sample_comparable_to_best$p_val_table,
            file = "tables/all_samples_comparable_to_best_pval.csv",
            row.names = FALSE)
      
  write.csv(x = all_sample_comparable_to_best$W_table,
            file = "tables/all_samples_comparable_to_best_W.csv",
            row.names = FALSE)
  
  all_sample_comparable_to_best$p_val_table %>%
    filter(pa_AUC > 0.05) -> all_value_useful_models
  
    
  full_output %>%
    filter(n_presence <= 20) %>%
    select(-n_presence,-n_background,-bg_method,-pres_method,-ratio_method,
           -n_pa_absence,-n_pa_presence,-species)%>%
    select(model,contains("pa_"))%>%
    sig_different_by_column(filter_by = "pa_AUC") -> ssss_sample_comparable_to_best
    
  
  full_output %>%
    filter(n_presence <= 20) %>%
    select(-n_presence,-n_background,-bg_method,-pres_method,-ratio_method,
           -n_pa_absence,-n_pa_presence,-species)%>%
    select(model,contains("pa_"))%>%
    sapply(function(y){ sum(length(which(is.na(y))))/length(y)})
  
    
  full_output %>%
    filter(n_presence <= 20) %>%
    select(-n_presence,-n_background,-bg_method,-pres_method,-ratio_method,
           -n_pa_absence,-n_pa_presence,-species)%>%
    select(model,contains("pa_"))%>%
    pivot_longer(cols = 2:10)
  
  sum(is.na(full_output$full_pAUC_specificity))/nrow(full_output)

  
  ssss_sample_comparable_to_best$p_val_table %>%
    filter(pa_AUC > 0.05)%>% mutate(model2=model)-> ssss_useful_models
  
      
##################
library(svglite)
  
  # PLotting overall models that are comparable to maxnet

      #Overall comparison
      
      full_output %>%
        filter(model != "CVmaxnet" ) %>%
        mutate(model = factor(model,
                              levels = full_output %>%
                                group_by(model) %>%
                                summarize(sm = mean(pa_sensitivity,na.rm = TRUE))%>%
                                arrange(sm) %>%
                                pull(model))) %>%
        filter(model %in% all_sample_comparable_to_maxnet ) %>%
        select(model,pa_AUC,pa_specificity,pa_sensitivity) %>%
        pivot_longer(contains("pa_"),names_to = "metric") %>%
        mutate(metric = gsub(pattern = "pa_",replacement = "",x =metric))%>%
        ggplot(mapping = aes(x=model,y=value))+
        geom_boxplot(notch = TRUE)+ #notch conveys additional info on median
        #geom_violin()+
        facet_wrap(~ metric,
                   ncol = 1,
                   scales = "free_y")+
        theme_bw()+
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              strip.text.x = element_text(size = 12))+
        xlab(NULL)+
        ylab(NULL)-> maxent_comparable_models_plot
      
      maxent_comparable_models_plot
      
      ggsave(plot = maxent_comparable_models_plot,
             filename = "figures/maxent_comparable_models_plot.svg",
             width = 10,height = 6,units = "in",dpi = 600)
      
      ggsave(plot = maxent_comparable_models_plot,
             filename = "figures/maxent_comparable_models_plot.jpg",
             width = 10,height = 6,units = "in",dpi = 600)
      
#################
  
# How many species with 20 or fewer records?
      
full_output %>%
        select(species,n_presence)%>%
        unique()%>%
        filter(n_presence <= 20) %>%
        summarise(n=n())
      
# How many species with more than 100 records?
      
      full_output %>%
        select(species,n_presence)%>%
        unique()%>%
        filter(n_presence > 100) %>%
        summarise(n=n())
      
      85*3*10
      
##########################################      

  # Model performance vs model type    

      
      full_output %>%
        filter(!is.na(bg_method)) %>%
        select(model,
               pres_method,
               bg_method,
               pa_AUC,
               pa_sensitivity,
               pa_specificity)%>%
        pivot_longer(cols = contains("pa_"),
                     names_to = "metric") %>%
        mutate(metric = gsub(pattern = "pa_",replacement = "",x=metric)) %>%
        # group_by(model,pres_method,bg_method,metric) %>%
        # summarise(mean_value = mean(value,na.rm = TRUE),
        #           ci_low = quantile(x=value, probs = 0.25,na.rm = TRUE),
        #           ci_high = quantile(x=value, probs = 0.75,na.rm = TRUE)) %>%
        ggplot(mapping = aes(y=value, x=model))+
        geom_hline(data = full_output %>%
                     filter(model == "maxnet") %>%
                     select(model,
                            pa_AUC,
                            pa_sensitivity,
                            pa_specificity)%>%
                     pivot_longer(cols = contains("pa_"),
                                  names_to = "metric") %>%
                     mutate(metric = gsub(pattern = "pa_",replacement = "",x=metric))%>%
                     group_by(model,metric)%>%
                     summarise(median = median(value,na.rm = TRUE)),
                   mapping = aes(yintercept = median),
                   lty = 2,color = "black"
        )+
        geom_boxplot(mapping = aes(color = pres_method),
                     notch = TRUE,
                     na.rm = TRUE,
                     fill="transparent",show.legend = FALSE)+
        facet_wrap(~ metric,ncol = 1
                   #,scales = "free_y"
                   )+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 6.5)),
          colour = "grey"
        )+
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 12.5)),
          colour = "grey"
        )+
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 18.5)),
          colour = "grey"
        )+
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 24.5)),
          colour = "grey"
        )+
        xlab(NULL)+
        ylab(NULL) -> pb_performance

      pb_performance
              
      ggsave(plot = pb_performance,
             filename = "figures/pb_performance_plot.svg",
             width = 10,
             height = 6,
             units = "in",
             dpi = 600)
      
      ggsave(plot = pb_performance,
             filename = "figures/pb_performance_plot.jpg",
             width = 10,
             height = 6,
             units = "in",
             dpi = 600)
      
###################################
      
    # Table of best performance by model type
      
  full_output %>%
        filter(model != "CVmaxnet" ) %>%
        group_by(species) %>%
        mutate(max_pa_AUC = max(pa_AUC,na.rm = TRUE)) %>%
        ungroup() %>%
        select(model,species,pa_AUC,max_pa_AUC) %>%
        arrange(species) %>%
        filter(max_pa_AUC == pa_AUC) %>%
        unique() %>%
        group_by(model) %>%
        summarise(times_won = n()) %>%
        mutate(pct_won = times_won/length(unique(full_output$species))*100)%>%
        arrange(-times_won)->times_won_table

      sum(times_won_table$pct_won)
      
# Table of best performance by model type, small sample size
      
      
      full_output %>%
        filter(model != "CVmaxnet" ) %>%
        filter(n_presence <= 20) %>%
        group_by(species) %>%
        mutate(max_pa_AUC = max(pa_AUC,na.rm = TRUE)) %>%
        ungroup() %>%
        select(model,species,pa_AUC,max_pa_AUC) %>%
        arrange(species) %>%
        filter(max_pa_AUC == pa_AUC) %>%
        unique() %>%
        mutate(n_focal_species = n())%>%
        group_by(model) %>%
        reframe(times_won = n(),
                  pct_won = n()/n_focal_species*100) %>%
        unique()%>%
        arrange(-times_won)->times_won_table_ssss

      sum(times_won_table_ssss$pct_won)
      
      full_output %>%
        filter(model != "CVmaxnet" ) %>%
        filter(n_presence <= 20)%>%
        pull(species)%>%
        unique()%>%
        length()

#####################################
    
    # Table of model poorness

poor_ssss_model_table <-      
full_output %>%
        filter(n_presence <= 20) %>%
        group_by(model)%>%
        summarise(pct_na_sens = round((sum(is.na(pa_sensitivity))/n())*100),
                  pct_na_spec = round((sum(is.na(pa_specificity))/n())*100),
                  pct_na_auc =round(( sum(is.na(pa_AUC))/n())*100),
                  pct_bad_aucs = round(((sum(pa_AUC < 0.5,na.rm = TRUE)+sum(is.na(pa_AUC)))/n())*100))

poor_ssss_model_table %>%
  arrange(pct_bad_aucs)%>%      
write.csv(file = "tables/poor_model_performance_table_ssss.csv",
          row.names = FALSE)      
      
poor_model_table <-      
  full_output %>%
#  filter(n_presence <= 20) %>%
  group_by(model)%>%
  summarise(pct_na_sens = round((sum(is.na(pa_sensitivity))/n())*100),
            pct_na_spec = round((sum(is.na(pa_specificity))/n())*100),
            pct_na_auc =round(( sum(is.na(pa_AUC))/n())*100),
            pct_bad_aucs = round(((sum(pa_AUC < 0.5,na.rm = TRUE)+sum(is.na(pa_AUC)))/n())*100))

############################

# Model runtime table

full_output %>%
  group_by(model) %>%
  summarise(median_runtime_s = median(runtime,na.rm = TRUE))%>%
  mutate(median_runtime_s = round(median_runtime_s,digits = 3)) %>%
  arrange(median_runtime_s) %>%
  write.csv("tables/median_model_runtimes.csv",row.names = FALSE)
