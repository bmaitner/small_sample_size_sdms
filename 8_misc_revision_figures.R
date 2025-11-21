# This file contains additional figures that were added during revision

#########################################

# Example map of sites
library(sf)
library(tidyverse)
library(disdat)

regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
data_i <- disData(region = "AWT")

data_i$bg %>%
  select(x,y) %>%
  unique()

data_i$po %>%
  select(x,y) %>%
  unique()


# Check whether any of the P/A sites are in the BG or PO

  any(data_i$pa$siteid %in% data_i$bg$siteid)
  any(data_i$pa$siteid %in% data_i$po$siteid)
  any(data_i$po$siteid %in% data_i$bg$siteid)
  
# Plot showing spatial distribution of one location

awt_plot <-   
data_i$env %>%
  select(x,y)%>%
  mutate(`Coordinate type` = "Presence/Absence") %>%
    st_as_sf(coords=c("x","y")) %>%
  bind_rows(data_i$bg %>%
              select(x,y)%>%
              mutate(`Coordinate type` = "Background") %>%
              st_as_sf(coords=c("x","y"))
            ) %>%
  bind_rows(data_i$po %>%
              select(x,y)%>%
              mutate(`Coordinate type` = "Presence Only") %>%
              st_as_sf(coords=c("x","y"))
            )%>%
  ggplot(mapping = aes(color=`Coordinate type`))+
  geom_sf()+
  geom_sf(data = . %>%
            filter(`Coordinate type` == "Presence/Absence"))+
  theme_void()+
  viridis::scale_color_viridis(discrete = TRUE)

ggsave(plot = awt_plot,
       filename = "figures/AWT_spatial_distribution_plot.jpg",
       height = 5,width = 4,units = "in",dpi = 300)



############################################################


# Plot showing consistency of sensitivity/specificity across settings

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


pnp_components <- c("rangebagging",
                    "gaussian",
                    "kde",
                    "vine",
                    "lobagoc")


models_to_evaluate <-
        data.frame(presence_method = pnp_components,
                   background_method = "none")

dr_models_to_evaluate <- c("ulsif","rulsif","maxnet")
                    

quantiles_to_evaluate <- c(0.05, 0.10, 0.15, 0.20, 0.25)


quantile_tempfile_full <- "outputs/temp_quantile_bakeoff_output_full.rds"



if(file.exists(quantile_tempfile_full)){
    
  quantile_variation_output <- readRDS(quantile_tempfile_full)
  
}else{
  
  quantile_variation_output <- NULL
  
}



for(q in 1:length(quantiles_to_evaluate)){
  
  quantile_q <- quantiles_to_evaluate[q]
  
  for(i in 1:nrow(models_to_evaluate)){

    if(any(quantile_variation_output$quantile == quantile_q &
           quantile_variation_output$model == paste(models_to_evaluate$presence_method[i],
                                                    models_to_evaluate$background_method[i],
                                                    sep = "/"))
       ){
      
      message("Skipping quantile ",quantile_q, " model ",
              paste(models_to_evaluate$presence_method[i],
                    models_to_evaluate$background_method[i],
                    sep = "/")," already done");
      
      next
      }
    
    out_i <- evaluate_disdat(presence_method = models_to_evaluate$presence_method[i],
                    background_method = models_to_evaluate$background_method[i],
                    verbose = TRUE,
                    ratio_method = NULL,
                    quantile = quantile_q,
                    ncl = 5)
    
    # pull out the relevant info (model, quantile, sens, spec, AUC)
    
                
      formatted_out_i <-
              out_i$full_model_stats %>%
                    select(species,
                           pa_sensitivity,
                           pa_specificity,
                           pa_AUC)%>%
                    mutate(model = paste(models_to_evaluate$presence_method[i],
                                         models_to_evaluate$background_method[i],
                                         sep = "/"),
                           quantile = quantile_q) %>%
        mutate(pa_AUC = as.numeric(pa_AUC))

      quantile_variation_output <-
        bind_rows(quantile_variation_output,
                  formatted_out_i)
                  
    # save temporary bits
    quantile_variation_output %>%
      saveRDS(file = file.path(quantile_tempfile_full))
    
    rm(out_i,formatted_out_i)
    
  }#pnp
  
  for(j in 1:length(dr_models_to_evaluate)){
    
    if(any(quantile_variation_output$quantile == quantile_q &
           quantile_variation_output$model == dr_models_to_evaluate[j])
    ){next}
    
    out_j <- evaluate_disdat(presence_method = NULL,
                             background_method = NULL,
                             verbose = TRUE,
                             ratio_method = dr_models_to_evaluate[j],
                             quantile = quantile_q,
                             ncl = 5)
    
    # pull out the relevant info (model, quantile, sens, spec, AUC)
    
    
    quantile_variation_output <-
      bind_rows(quantile_variation_output,
                out_j$full_model_stats %>%
                  select(species,
                         pa_sensitivity,
                         pa_specificity,
                         pa_AUC)%>%
                  mutate(model = dr_models_to_evaluate[j],
                         quantile = quantile_q) %>%
                  mutate(pa_AUC = as.numeric(pa_AUC)))
    
    
    # save temporary bits
    
    quantile_variation_output %>%
      saveRDS(file = file.path(quantile_tempfile_full))
    
    rm(out_j,formatted_out_j)
    
    
  }#dr
  
  
  # save temporary bits
  
    quantile_variation_output %>%
      saveRDS(file = file.path(quantile_tempfile_full))

} #q loop


library(ggrepel)

# Plot
  # sens vs spec, color by model, facet by quantile
  # alternatively, can calculate rank, then plot median rank +/- max/min



presence_counts <- 
  readRDS("outputs/full_model_output_all.RDS") %>%
  select(species,
         n_presence)%>%
  na.omit()%>%
  unique()

quantile_variation_output <- quantile_variation_output %>%
  left_join(presence_counts,relationship = "many-to-many")


quantile_variation_output %>%
  #filter(n_presence <= 20)%>%
  # filter(method %in% c("kde / none","vine / none","gaussian / none",
  #                      "maxnet","ulsif","rulsif",
  #                      "rangebagging / none","lobagoc / none"))%>%
  mutate(model = gsub(pattern = "/none",replacement = "",x=model))%>%
  mutate(model = gsub(pattern = "maxnet",replacement = "Maxnet",x=model))%>%
  mutate(model = gsub(pattern = "vine",replacement = "Vine",x=model))%>%
  mutate(model = gsub(pattern = "gaussian",replacement = "Gaussian",x=model))%>%
  mutate(model = gsub(pattern = "kde",replacement = "KDE",x=model))%>%
  mutate(model = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=model))%>%
  mutate(model = gsub(pattern = "ulsif",replacement = "uLSIF",x=model))%>%
  mutate(model = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=model)) %>%  
  #mutate(sens_spec_ratio = pa_specificity-pa_sensitivity) %>%
  group_by(model,quantile)%>%
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
                       label = model,
                       color = model))+
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
                  point.padding = 0)+
  facet_wrap(~quantile) -> svs_variation

svs_variation


  
quantile_variation_rankings <-
quantile_variation_output %>%
  #filter(n_presence <= 20)%>%
  # filter(method %in% c("kde / none","vine / none","gaussian / none",
  #                      "maxnet","ulsif","rulsif",
  #                      "rangebagging / none","lobagoc / none"))%>%
  mutate(model = gsub(pattern = "/none",replacement = "",x=model))%>%
  mutate(model = gsub(pattern = "maxnet",replacement = "Maxnet",x=model))%>%
  mutate(model = gsub(pattern = "vine",replacement = "Vine",x=model))%>%
  mutate(model = gsub(pattern = "gaussian",replacement = "Gaussian",x=model))%>%
  mutate(model = gsub(pattern = "kde",replacement = "KDE",x=model))%>%
  mutate(model = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=model))%>%
  mutate(model = gsub(pattern = "ulsif",replacement = "uLSIF",x=model))%>%
  mutate(model = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=model)) %>%  
  #mutate(sens_spec_ratio = pa_specificity-pa_sensitivity) %>%
  group_by(model,quantile)%>%
  #summarise(sens_spec_ratio = median(na.omit(sens_spec_ratio)))%>%
  summarise(mean_sensitivity = mean(na.omit(pa_sensitivity)),
            mean_specificity = mean(na.omit(pa_specificity)),
            median_sensitivity = median(na.omit(pa_sensitivity)),
            median_specificity = median(na.omit(pa_specificity)),
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
  arrange(-mean_sensitivity)%>%
    group_by(quantile)%>%
  mutate(sens_rank = row_number())%>%
    arrange(-mean_specificity)%>%
    group_by(quantile)%>%
  mutate(spec_rank = row_number())%>%
  group_by(quantile)


sens_to_thresholding <-
ggplot(data = quantile_variation_rankings,
       aes(x = quantile,
           y = mean_specificity/mean_sensitivity,
           color = model)) +
  geom_point()+
  geom_line()+
  scale_y_log10()+
  theme_bw()+
  ylab("Mean Specificity / Mean Sensitivity")


ggsave(plot = sens_to_thresholding,
       filename = "figures/thresholding_impact_on_sensitivity_specificity.jpg",
       width = 4,
       height = 3,
       units = "in",
       dpi = 300)


################################################################################

  # Comparison with Valavi et al 2021

library(ggplot2)
library(tidyverse)
library(ggrepel)

# Get Valavi data

  source("R/get_valavi_stats.R")
  
  if(!file.exists("data/manual_downloads/Valavi/Valavi_model_stats.RDS")){
  
    valavi_stats <- get_valavi_stats(models_prediction_folder = "data/manual_downloads/Valavi/Models_prediction/")
    
    valavi_stats %>% saveRDS(file = "data/manual_downloads/Valavi/Valavi_model_stats.RDS")
      
  }else{
  
    valavi_stats <- readRDS(file = "data/manual_downloads/Valavi/Valavi_model_stats.RDS")
    
  }

# Load our data

combined_stats <- 
  readRDS("outputs/full_model_output_all.RDS") %>%
    mutate(pa_AUC = as.numeric(pa_AUC))%>%
  dplyr::select(spid = species,
                model = method,
                pa_AUC,
                pa_correlation) %>%
    mutate(author = "Maitner et al.") %>%
    bind_rows(valavi_stats %>%
                mutate( author = "Valavi et al.")%>%
                select(-region))

presence_counts <- 
  readRDS("outputs/full_model_output_all.RDS") %>%
  select(spid = species,
         n_presence)%>%
  na.omit()%>%
  unique()

combined_stats <- combined_stats %>%
  left_join(presence_counts,relationship = "many-to-many")

# Ranking by AUC

combined_stats_data_poor_summary <-
  combined_stats %>%
  filter(n_presence <= 20)%>%
  #mutate(model = gsub(pattern = " / none",replacement = "",x=model))%>%
  mutate(model = gsub(pattern = "maxnet",replacement = "Maxnet",x=model))%>%
  mutate(model = gsub(pattern = "vine",replacement = "Vine",x=model))%>%
  mutate(model = gsub(pattern = "gaussian",replacement = "Gaussian",x=model))%>%
  mutate(model = gsub(pattern = "kde",replacement = "KDE",x=model))%>%
  mutate(model = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=model))%>%
  mutate(model = gsub(pattern = "ulsif",replacement = "uLSIF",x=model))%>%
  mutate(model = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=model))%>%
  group_by(author,model)%>%
  summarise(median_AUC = median(na.omit(pa_AUC)),
            mean_AUC = mean(pa_AUC,na.rm=TRUE),
            SD_AUC = sd(pa_AUC,na.rm=TRUE),
            median_Corr = median(na.omit(pa_correlation)),
            mean_Corr = mean(na.omit(pa_correlation)),
            SD_Corr = sd(na.omit(pa_correlation)),
            AUC_Q1 = quantile(x = pa_AUC,
                              probs = 0.25,
                              na.rm=TRUE),
            AUC_Q3 = quantile(x = pa_AUC,
                              probs = 0.75,
                              na.rm=TRUE),
            correlation_Q1 = quantile(x = pa_correlation,
                                      probs = 0.25,
                                      na.rm=TRUE),
            correlation_Q3 = quantile(x = pa_correlation,
                                      probs = 0.75,
                                      na.rm=TRUE)
            
  ) %>%
  ungroup() %>%
  arrange(-median_AUC) |>
  mutate(across(median_AUC:correlation_Q3, ~round(.x,digits = 3))) |>
  pivot_longer(
    cols = matches("^(mean|SD)_"),
    names_to = c(".value", "var"),
    names_pattern = "(mean|SD)_(.*)"
  ) %>%
  mutate(mean = sprintf("%s (%s)", mean, SD)) %>%
  select(-SD) %>%
  pivot_wider(
    names_from = var,
    values_from = mean,
    names_prefix = "mean "
  )|>
  relocate(author,model,median_AUC,`mean AUC`,`mean Corr`)

  write.csv(x =   combined_stats_data_poor_summary,
            file = "tables/small_sample_size_comparison_to_Valavi_et_al.csv",
            row.names = FALSE)


  combined_stats_all_spp_summary <-
    combined_stats %>%
    #mutate(model = gsub(pattern = " / none",replacement = "",x=model))%>%
    mutate(model = gsub(pattern = "maxnet",replacement = "Maxnet",x=model))%>%
    mutate(model = gsub(pattern = "vine",replacement = "Vine",x=model))%>%
    mutate(model = gsub(pattern = "gaussian",replacement = "Gaussian",x=model))%>%
    mutate(model = gsub(pattern = "kde",replacement = "KDE",x=model))%>%
    mutate(model = gsub(pattern = "rangebagging",replacement = "Range-Bagging",x=model))%>%
    mutate(model = gsub(pattern = "ulsif",replacement = "uLSIF",x=model))%>%
    mutate(model = gsub(pattern = "lobagoc",replacement = "LOBAG-OC",x=model))%>%
    group_by(author,model)%>%
    summarise(median_AUC = median(na.omit(pa_AUC)),
              mean_AUC = mean(pa_AUC,na.rm=TRUE),
              SD_AUC = sd(pa_AUC,na.rm=TRUE),
              median_Corr = median(na.omit(pa_correlation)),
              mean_Corr = mean(na.omit(pa_correlation)),
              SD_Corr = sd(na.omit(pa_correlation)),
              AUC_Q1 = quantile(x = pa_AUC,
                                probs = 0.25,
                                na.rm=TRUE),
              AUC_Q3 = quantile(x = pa_AUC,
                                probs = 0.75,
                                na.rm=TRUE),
              correlation_Q1 = quantile(x = pa_correlation,
                                        probs = 0.25,
                                        na.rm=TRUE),
              correlation_Q3 = quantile(x = pa_correlation,
                                        probs = 0.75,
                                        na.rm=TRUE)
              
    ) %>%
    ungroup() %>%
    arrange(-median_AUC) |>
    mutate(across(median_AUC:correlation_Q3, ~round(.x,digits = 3))) |>
    pivot_longer(
      cols = matches("^(mean|SD)_"),
      names_to = c(".value", "var"),
      names_pattern = "(mean|SD)_(.*)"
    ) %>%
    mutate(mean = sprintf("%s (%s)", mean, SD)) %>%
    select(-SD) %>%
    pivot_wider(
      names_from = var,
      values_from = mean,
      names_prefix = "mean "
    )|>
    relocate(author,model,median_AUC,`mean AUC`,`mean Corr`)
  
  write.csv(x =   combined_stats_all_spp_summary,
            file = "tables/all_sample_size_comparison_to_Valavi_et_al.csv",
            row.names = FALSE)
  
  

# Ranking by number of "wins"
  
  auc_wins_data_poor_summary <-
  combined_stats %>%
    filter(n_presence <= 20)%>%
    group_by(spid) %>%
    arrange(desc(pa_AUC)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(author,model)%>%
    summarize(n_wins = n())%>%
    arrange(-n_wins)

  auc_wins_all_summary <-
  combined_stats %>%
    group_by(spid) %>%
    arrange(desc(pa_AUC)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(author,model)%>%
    summarize(n_wins = n())%>%
    ungroup()%>%
    mutate(pct_wins = (n_wins/sum(n_wins))*100)%>%
    arrange(-n_wins)
  
  # How often does one of our models win?
  
  combined_stats %>%
    #filter(model != "maxnet") %>%
    group_by(spid) %>%
    arrange(desc(pa_AUC)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(author)%>%
    summarize(n_wins = n())%>%
    arrange(-n_wins)%>% #93 times (vs 133).  (88 if you exclude maxnet)
    ungroup() %>%
    mutate(pct_wins = n_wins/sum(n_wins)*100)

  combined_stats %>%
    filter(n_presence <= 20)%>%
    group_by(spid) %>%
    arrange(desc(pa_AUC)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(author)%>%
    summarize(n_wins = n())%>%
    arrange(-n_wins) %>% # 19 times (vs 15)
    ungroup() %>%
    mutate(pct_wins = n_wins/sum(n_wins)*100)
  
  
  
  

# Testing which models have different distributions

  # small sample sizess only
  
  for(i in 1:length(unique(combined_stats$model))){
    for(j in 1:length(unique(combined_stats$model))){
      
      model_i <- unique(combined_stats$model)[i]
      model_j <- unique(combined_stats$model)[j]
      
      data_x <- combined_stats %>%
        filter(n_presence <= 20) %>%
        filter(model==model_i)%>%
        pull(pa_AUC)
      
      data_y <- combined_stats %>%
        filter(n_presence <= 20) %>%
        filter(model==model_j)%>%
        pull(pa_AUC)
      
      ti <-tryCatch(wilcox.test(x = data_x,
                                y = data_y),
                    error = function(e){e})
      
      if(inherits(ti,"error")){
        ti$p.value <- NA
        ti$statistic <- NA
        
      }
      
      out_i <- data.frame(model = model_i,
                          comparison = model_j,
                          pval = ti$p.value,
                          W = ti$statistic,
                          comparison_difference = mean(data_x,na.rm = TRUE) - mean(data_y,na.rm = TRUE))
      
      if(i==1 & j==1){ttest_small_out <- out_i}else(ttest_small_out <- bind_rows(ttest_small_out,out_i))
      
      rm(out_i,ti,model_i,model_j,data_x,data_y)

    }#j
  } #i 

  
ttest_small_matrix <-    
ttest_small_out %>%
  select(model,comparison,pval)%>%
  mutate(sig_dif = (pval <= 0.05) %>% as.numeric() )%>%
  select(model,comparison,sig_dif)%>%
  pivot_wider(values_from = sig_dif,names_from = comparison) %>%
  column_to_rownames("model")

write.csv(ttest_small_matrix,
          "tables/small_sample_size_model_comparison_sig_test_matrix.csv")

write.csv(ttest_small_out,
          "tables/small_sample_size_model_comparison.csv",
          row.names = FALSE)



# all sample sizes

for(i in 1:length(unique(combined_stats$model))){
  for(j in 1:length(unique(combined_stats$model))){
    
    model_i <- unique(combined_stats$model)[i]
    model_j <- unique(combined_stats$model)[j]
    
    data_x <- combined_stats %>%
      #filter(n_presence <= 20) %>%
      filter(model==model_i)%>%
      pull(pa_AUC)
    
    data_y <- combined_stats %>%
      #filter(n_presence <= 20) %>%
      filter(model==model_j)%>%
      pull(pa_AUC)
    
    ti <-tryCatch(wilcox.test(x = data_x,
                              y = data_y),
                  error = function(e){e})
    
    if(inherits(ti,"error")){
      ti$p.value <- NA
      ti$statistic <- NA
      
    }
    
    out_i <- data.frame(model = model_i,
                        comparison = model_j,
                        pval = ti$p.value,
                        W = ti$statistic,
                        comparison_difference = mean(data_x,na.rm = TRUE) - mean(data_y,na.rm = TRUE))
    
    if(i==1 & j==1){ttest_all_out <- out_i}else(ttest_all_out <- bind_rows(ttest_all_out,out_i))
    
    rm(out_i,ti,model_i,model_j,data_x,data_y)
    
  }#j
} #i 


ttest_all_matrix <-    
  ttest_all_out %>%
  select(model,comparison,pval)%>%
  mutate(sig_dif = (pval <= 0.05) %>% as.numeric() )%>%
  select(model,comparison,sig_dif)%>%
  pivot_wider(values_from = sig_dif,names_from = comparison) %>%
  column_to_rownames("model")

# How many models are garbage?

Valavi_comparison_bad_model_small_sample_size <-
combined_stats %>%
  filter(n_presence <= 20)%>%
  group_by(model) %>%
  mutate(bad_model = (pa_AUC <= 0.5) %>% as.numeric())%>%
  mutate(bad_model = case_when(is.na(bad_model) ~ 1,.default = bad_model))%>%
  group_by(author,model)%>%
  summarize(pct_bad_models = (sum(bad_model)/n())*100)%>%
  ungroup()%>%
  arrange(pct_bad_models)

write.csv(Valavi_comparison_bad_model_small_sample_size,
          "tables/Valavi_comparison_bad_model_small_sample_size.csv",
          row.names = FALSE)

# presence informaiton

min(presence_counts$n_presence)#2
median(presence_counts$n_presence)#57
mean(presence_counts$n_presence)#233
max(presence_counts$n_presence)#5822

###################################

# Investigate overlap between model predictions

  # Calculate overlap stats

    if(file.exists("outputs/model_prediction_agreement.RDS")){
      
      model_prediction_agreement <- readRDS("outputs/model_prediction_agreement.RDS")
      
    }else{
      
      library(arrow)
      source("R/get_prediction_overlap.R")
      
      
      
      model_prediction_agreement <-
        get_prediction_overlap(dataset_folder = "outputs/model_predictions/",
                               self_comparison = TRUE)
      
      
      saveRDS(object = model_prediction_agreement,
              file = "outputs/model_prediction_agreement.RDS")
      
      
      
    }

  # Calc model sens-spec distance

  combined_stats <- 
    readRDS("outputs/full_model_output_all.RDS") %>%
    mutate(pa_AUC = as.numeric(pa_AUC)) %>%
    mutate(method = gsub(pattern = " ",replacement = "",x = method))

  source("R/get_sensspec_distance.R")

    sens_spec_dist <- get_sensspec_distance(combined_stats = combined_stats,
                                            self_comparison = TRUE,
                                            model_pairs = model_prediction_agreement %>%
                                              ungroup%>%
                                              select(model_a,model_b)%>%
                                              unique())
    
    
    
    


  # Join sens-spec distance to agreement

    sensspec_v_agreement<-
    sens_spec_dist %>%
      full_join(model_prediction_agreement,
                by = c("species","model_a","model_b"))

    
      
  # Visualize

    library(ggplot2)
    library(ggpubr)
    
    sensspec_v_agreement %>%
      ggplot(mapping = aes(x=sens_spec_distance,
                           y = model_agreement)) +
      geom_point()
    
    
    sensspec_v_agreement %>%
      ggplot(mapping = aes(x = sens_spec_distance,
                           y = Jaccard)) +
      geom_point()
    
    
    sensspec_v_agreement %>%
      group_by(model_a,model_b)%>%
      summarise(sens_spec_distance = mean(sens_spec_distance,na.rm=TRUE),
                model_agreement = mean(model_agreement,na.rm=TRUE))%>%
      ggplot(mapping = aes(x=sens_spec_distance,
                           y = model_agreement)) +
      geom_point()

    
    sensspec_v_agreement %>%
      group_by(model_a,model_b)%>%
      summarise(sens_spec_distance = mean(sens_spec_distance,na.rm=TRUE),
                Jaccard = mean(Jaccard,na.rm=TRUE))%>%
      ggplot(mapping = aes(x=sens_spec_distance,
                           y = Jaccard)) +
      geom_point()
    
    sensspec_vagreement_jaccard_plot <-
      sensspec_v_agreement %>%
      group_by(model_a,model_b)%>%
      summarise(sens_spec_distance = mean(sens_spec_distance,na.rm=TRUE),
                Jaccard = mean(Jaccard,na.rm=TRUE))%>%
      ggplot(mapping = aes(x=sens_spec_distance,
                           y = Jaccard)) +
      geom_point()+
      stat_cor(method = "pearson",label.x = 0.7)+
      theme_bw()+
      xlab("Sensitivity/Specificity Distance")+
      ylab("Jaccard Similarity")
    
    sensspec_vagreement_plot <-
    sensspec_v_agreement %>%
      group_by(model_a,model_b)%>%
      summarise(sens_spec_distance = mean(sens_spec_distance,na.rm=TRUE),
                model_agreement = mean(model_agreement,na.rm=TRUE))%>%
      ggplot(mapping = aes(x=sens_spec_distance,
                           y = model_agreement)) +
      geom_point()+
      stat_cor(method = "pearson",label.x = 0.7)+
      theme_bw()+
      xlab("Sensitivity/Specificity Distance")+
      ylab("Model Agreement")
    
    ggsave(filename = "figures/sensspec_vagreement_plot.jpeg",
           plot = sensspec_vagreement_plot,
           dpi = 600,
           width = 5,
           height = 5,
           units = "in")
    
    ggsave(filename = "figures/sensspec_vjaccard_plot.jpeg",
           plot = sensspec_vagreement_jaccard_plot,
           dpi = 600,
           width = 5,
           height = 5,
           units = "in")

    combined_agreement_plot <-
    ggarrange(sensspec_vagreement_plot,
              sensspec_vagreement_jaccard_plot,
              ncol = 1,
              labels = c("A","B"))    
    

    ggsave(filename = "figures/sensspec_v_agreement_and_jaccard_plot.jpeg",
           plot = combined_agreement_plot,
           dpi = 600,
           width = 5,
           height = 10,
           units = "in")
    
##################################################

    # get counts of predictors used in each dataset    
    
    library(disdat)

    regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
    
    pred_counts <- NULL
    
    for (i in 1:length(regions)){
      
      
      region_i <- regions[i]
      data_i <- disData(region = region_i)
      
      message(paste("Starting region ",i, " of ", length(regions),": ",region_i))
      
      
      #Remove categorical predictors  
      
      if(region_i == "CAN"){
        
        data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("ontveg"))] 
        data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("ontveg"))]
        data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("ontveg"))]
        
      }
      
      if(region_i == "NSW"){
        data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("disturb","soilfert","vegsys"))] 
        data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("disturb","soilfert","vegsys"))]
        data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("disturb","soilfert","vegsys"))]
      }
      
      if(region_i == "NZ"){
        data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("age","toxicats"))] 
        data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("age","toxicats"))]
        data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("age","toxicats"))]
      }
      
      if(region_i == "SWI"){
        data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("calc","sfroyy"))] 
        data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("calc","sfroyy"))]
        data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("calc","sfroyy"))]
      }
      
      pred_counts <- bind_rows(pred_counts,
                       data.frame(region = region_i,
                                  n_predictors = data_i$bg %>%
                                    select(7:ncol(.)) %>%
                                    ncol()
                                  )
                       )
                       
    }#end i
    
    rm(regions)  

    
    pred_counts    
    
#############################################################    
# BIEN occurrence accuracy
    
  occ_counts <- readRDS("data/manual_downloads/BIEN_occs/occ_counts.RDS")    
  
  ssss_counts <- occ_counts %>%
    filter(n_occs <= 100)
  
  
  library(BIEN)  

    
    fl_species <- BIEN_list_state(country = "United States",
                                  state = "Florida")
    
    ssss_counts %>%
      mutate(species = gsub(pattern = "_",
                            replacement = " ",
                            x = species)) -> ssss_counts
    
    fl_counts <-
      ssss_counts %>%
      filter(species %in% fl_species$scrubbed_species_binomial)
    
    # Get coords 
    
    fl_occs <- BIEN:::.BIEN_sql(query = "SELECT scrubbed_species_binomial,
                     latitude,
                     longitude,
                     coord_uncertainty_m,
                     coordinate_inherent_uncertainty_m
                     FROM view_full_occurrence_individual
                     WHERE country in ( 'United States' )
                     AND state_province in ( 'Florida' )
                     AND scrubbed_species_binomial IS NOT NULL
                     AND (is_cultivated_observation = 0 OR is_cultivated_observation IS NULL)
                     AND is_location_cultivated IS NULL
                     AND (is_introduced=0 OR is_introduced IS NULL)
                     AND higher_plant_group NOT IN ('Algae','Bacteria','Fungi')
                     AND is_geovalid = 1
                     AND (georef_protocol is NULL OR georef_protocol<>'county centroid')
                     AND (is_centroid IS NULL OR is_centroid=0)
                     AND observation_type IN ('plot','specimen','literature','checklist')
                     AND scrubbed_species_binomial IS NOT NULL ;",
                     fetch.query = FALSE)
    
    fl_occs <- fl_occs %>%
      filter(scrubbed_species_binomial %in% fl_counts$species)
    
    
    fl_occs %>%
      summarise(mean_coord_uncertainty_m = mean(coord_uncertainty_m,na.rm=TRUE),
                mean_inherent_uncertainty_m = mean(coordinate_inherent_uncertainty_m,na.rm=TRUE),
                median_coord_uncertainty_m = median(coord_uncertainty_m,na.rm=TRUE),
                median_inherent_uncertainty_m = median(coordinate_inherent_uncertainty_m,na.rm=TRUE))
    
    # mean_coord_uncertainty_m mean_inherent_uncertainty_m median_coord_uncertainty_m
    #                  5861.967                    1010.975                        651
    # median_inherent_uncertainty_m
    #                      0.9767545 
    
    

    
    
      
###########################################################################

library(ggplot2)    
        
# Plots vs maxent    
    
    ga_vs_mn <-
      combined_stats %>%
      select(method,pres_method,species,pa_AUC,n_presence) %>%
      mutate(small_sample_size = case_when(n_presence <= 20 ~ "yes",
                                           n_presence > 20 ~ "no"))%>%
      filter(pres_method %in% "gaussian" |
               method %in% "maxnet") %>%
      left_join(x=.,
                y = dplyr::filter(.data = .,method %in% "maxnet")%>%
                  rename(method_comp=method,
                         pa_AUC_comp = pa_AUC)%>%
                  select(-pres_method,-small_sample_size),by="species",copy = TRUE) %>%
      filter(method != "maxnet")%>%
      na.omit()%>%
      ggplot(mapping = aes(x=pa_AUC_comp,y = pa_AUC,color=small_sample_size))+
      geom_point(alpha=0.5)+
      facet_wrap(~method)+
      xlab("Maxnet P/A AUC")+
      ylab("Focal P/A AUC")+
      geom_hline(yintercept = .5,linetype = 2)+
      geom_vline(xintercept = .5,linetype = 2)+
      geom_smooth(method = "lm",linewidth = 2)+
      geom_abline(slope = 1,intercept = 0)+
      labs(color="Small Sample Size")+
      theme_bw()
    
    kde_vs_mn <-
      combined_stats %>%
      select(method,pres_method,species,pa_AUC,n_presence) %>%
      mutate(small_sample_size = case_when(n_presence <= 20 ~ "yes",
                                           n_presence > 20 ~ "no"))%>%
      filter(pres_method %in% "kde" |
               method %in% "maxnet") %>%
      left_join(x=.,
                y = dplyr::filter(.data = .,method %in% "maxnet")%>%
                  rename(method_comp=method,
                         pa_AUC_comp = pa_AUC)%>%
                  select(-pres_method,-small_sample_size),by="species",copy = TRUE) %>%
      filter(method != "maxnet")%>%
      na.omit()%>%
      ggplot(mapping = aes(x=pa_AUC_comp,y = pa_AUC,color=small_sample_size))+
      geom_point(alpha=0.5)+
      facet_wrap(~method)+
      xlab("Maxnet P/A AUC")+
      ylab("Focal P/A AUC")+
      geom_hline(yintercept = .5,linetype = 2)+
      geom_vline(xintercept = .5,linetype = 2)+
      geom_smooth(method = "lm",linewidth = 2)+
      geom_abline(slope = 1,intercept = 0)+
      labs(color="Small Sample Size")+
      theme_bw()
    
    lb_vs_mn <-
      combined_stats %>%
      select(method,pres_method,species,pa_AUC,n_presence) %>%
      mutate(small_sample_size = case_when(n_presence <= 20 ~ "yes",
                                           n_presence > 20 ~ "no"))%>%
      filter(pres_method %in% "lobagoc" |
               method %in% "maxnet") %>%
      left_join(x=.,
                y = dplyr::filter(.data = .,method %in% "maxnet")%>%
                  rename(method_comp=method,
                         pa_AUC_comp = pa_AUC)%>%
                  select(-pres_method,-small_sample_size),by="species",copy = TRUE) %>%
      filter(method != "maxnet")%>%
      na.omit()%>%
      ggplot(mapping = aes(x=pa_AUC_comp,y = pa_AUC,color=small_sample_size))+
      geom_point(alpha=0.5)+
      facet_wrap(~method)+
      xlab("Maxnet P/A AUC")+
      ylab("Focal P/A AUC")+
      geom_hline(yintercept = .5,linetype = 2)+
      geom_vline(xintercept = .5,linetype = 2)+
      geom_smooth(method = "lm",linewidth = 2)+
      geom_abline(slope = 1,intercept = 0)+
      labs(color="Small Sample Size")+
      theme_bw()
    
    
    rb_vs_mn <-
      combined_stats %>%
      select(method,pres_method,species,pa_AUC,n_presence) %>%
      mutate(small_sample_size = case_when(n_presence <= 20 ~ "yes",
                                           n_presence > 20 ~ "no"))%>%
      filter(pres_method %in% "rangebagging" |
               method %in% "maxnet") %>%
      left_join(x=.,
                y = dplyr::filter(.data = .,method %in% "maxnet")%>%
                  rename(method_comp=method,
                         pa_AUC_comp = pa_AUC)%>%
                  select(-pres_method,-small_sample_size),by="species",copy = TRUE) %>%
      filter(method != "maxnet")%>%
      na.omit()%>%
      ggplot(mapping = aes(x=pa_AUC_comp,y = pa_AUC,color=small_sample_size))+
      geom_point(alpha=0.5)+
      facet_wrap(~method)+
      xlab("Maxnet P/A AUC")+
      ylab("Focal P/A AUC")+
      geom_hline(yintercept = .5,linetype = 2)+
      geom_vline(xintercept = .5,linetype = 2)+
      geom_smooth(method = "lm",linewidth = 2)+
      geom_abline(slope = 1,intercept = 0)+
      labs(color="Small Sample Size")+
      theme_bw()
    
    
    vn_vs_mn <-
      combined_stats %>%
      select(method,pres_method,species,pa_AUC,n_presence) %>%
      mutate(small_sample_size = case_when(n_presence <= 20 ~ "yes",
                                           n_presence > 20 ~ "no"))%>%
      filter(pres_method %in% "vine" |
               method %in% "maxnet") %>%
      left_join(x=.,
                y = dplyr::filter(.data = .,method %in% "maxnet")%>%
                  rename(method_comp=method,
                         pa_AUC_comp = pa_AUC)%>%
                  select(-pres_method,-small_sample_size),by="species",copy = TRUE) %>%
      filter(method != "maxnet")%>%
      na.omit()%>%
      ggplot(mapping = aes(x=pa_AUC_comp,y = pa_AUC,color=small_sample_size))+
      geom_point(alpha=0.5)+
      facet_wrap(~method)+
      xlab("Maxnet P/A AUC")+
      ylab("Focal P/A AUC")+
      geom_hline(yintercept = .5,linetype = 2)+
      geom_vline(xintercept = .5,linetype = 2)+
      geom_smooth(method = "lm",linewidth = 2)+
      geom_abline(slope = 1,intercept = 0)+
      labs(color="Small Sample Size")+
      theme_bw()
    
    
    dr_vs_mn <-
      combined_stats %>%
      select(method,ratio_method,species,pa_AUC,n_presence) %>%
      mutate(small_sample_size = case_when(n_presence <= 20 ~ "yes",
                                           n_presence > 20 ~ "no"))%>%
      filter(!is.na(ratio_method)) %>%
      left_join(x=.,
                y = dplyr::filter(.data = .,method %in% "maxnet")%>%
                  rename(method_comp=method,
                         pa_AUC_comp = pa_AUC)%>%
                  select(-ratio_method,-small_sample_size),by="species",copy = TRUE) %>%
      filter(method != "maxnet")%>%
      na.omit()%>%
      ggplot(mapping = aes(x=pa_AUC_comp,y = pa_AUC,color=small_sample_size))+
      geom_point(alpha=0.5)+
      facet_wrap(~method)+
      xlab("Maxnet P/A AUC")+
      ylab("Focal P/A AUC")+
      geom_hline(yintercept = .5,linetype = 2)+
      geom_vline(xintercept = .5,linetype = 2)+
      geom_smooth(method = "lm",linewidth = 2)+
      geom_abline(slope = 1,intercept = 0)+
      labs(color="Small Sample Size")+
      theme_bw()
    
    ggsave(plot = ga_vs_mn,
           filename = "figures/maxnet_v_gaussian.jpg",
           height = 6,width = 8,units = "in",dpi = 300)
    
    ggsave(plot = kde_vs_mn,
           filename = "figures/maxnet_v_kde.jpg",
           height = 6,width = 8,units = "in",dpi = 300)
    
    ggsave(plot = lb_vs_mn,
           filename = "figures/maxnet_v_lobagoc.jpg",
           height = 6,width = 8,units = "in",dpi = 300)
    
    ggsave(plot = rb_vs_mn,
           filename = "figures/maxnet_v_rangebagging.jpg",
           height = 6,width = 8,units = "in",dpi = 300)
    
    ggsave(plot = vn_vs_mn,
           filename = "figures/maxnet_v_vine.jpg",
           height = 6,width = 8,units = "in",dpi = 300)
    
    ggsave(plot = dr_vs_mn,
           filename = "figures/maxnet_v_density_ratio.jpg",
           height = 6,width = 8,units = "in",dpi = 300)
    
##############################################
    
  # Distribution of occurrence records for small sample size species in Elith dataset  
    
    combined_stats %>%
        filter(method=="maxnet")%>%
        filter(n_presence < 20) %>%
        pull(n_presence) %>% hist(breaks = 20) 
    
###############################################
    
  # AUC mean vs variance (after https://onlinelibrary.wiley.com/doi/10.1111/j.1472-4642.2008.00482.x )
    
library(ggrepel)    

    
small_sample_size_med_v_var <-        
combined_stats %>%
      filter(n_presence <= 20) %>%
      mutate(model_group = case_when(!is.na(pres_method) ~ pres_method,
                                     is.na(pres_method) ~ ratio_method))%>%
      group_by(method,model_group) %>%
      summarise(median_pa_auc = median(pa_AUC,na.rm = TRUE),
                variance_pa_auc = var(pa_AUC,na.rm = TRUE)) %>%
      ggplot(mapping = aes(x=median_pa_auc,y=variance_pa_auc,color=model_group))+
      geom_point()+
  geom_text_repel(aes(label = method),
                  min.segment.length = .1,
                  show.legend = FALSE,
                  max.overlaps=100)+
  theme_bw()+
  labs(color = "Algorithm Group")+
  xlab("Median P/A AUC")+
  ylab("Variance P/A AUC")
  #+scale_x_continuous(expand = c(0.1,0.1)) +
  #scale_y_continuous(expand = c(0.005,0.005))

ggsave(plot = small_sample_size_med_v_var,
       filename = "figures/auc_med_v_var_small_sample_size.jpg",
       width = 12,
       height = 10,
       units = "in",
       dpi = 300)


all_sample_size_med_v_var <-        
  combined_stats %>%
#  filter(n_presence <= 20) %>%
  mutate(model_group = case_when(!is.na(pres_method) ~ pres_method,
                                 is.na(pres_method) ~ ratio_method))%>%
  group_by(method,model_group) %>%
  summarise(median_pa_auc = median(pa_AUC,na.rm = TRUE),
            variance_pa_auc = var(pa_AUC,na.rm = TRUE)) %>%
  ggplot(mapping = aes(x=median_pa_auc,y=variance_pa_auc,color=model_group))+
  geom_point()+
  geom_text_repel(aes(label = method),
                  min.segment.length = .1,
                  show.legend = FALSE,
                  max.overlaps=100)+
  theme_bw()+
  labs(color = "Algorithm Group")+
  xlab("Median P/A AUC")+
  ylab("Variance P/A AUC")
#+scale_x_continuous(expand = c(0.1,0.1)) +
#scale_y_continuous(expand = c(0.005,0.005))


ggsave(plot = all_sample_size_med_v_var,
       filename = "figures/auc_med_v_var_all_sample_size.jpg",
       width = 12,
       height = 10,
       units = "in",
       dpi = 300)

