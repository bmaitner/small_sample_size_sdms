# 3a pnp bakeoff analyses

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

# Load and format cv data
  
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
      group_by(model) %>%
      summarise('median PA AUC' = na.omit(pa_AUC) %>% median(),
                'mean PA AUC' = na.omit(pa_AUC) %>% mean(),
                'mean PA sensitivity' = na.omit(pa_sensitivity) %>% mean(),
                'mean PA specificity' = na.omit(pa_specificity) %>% mean(),
                'mean PA pAUC (sensitivity 0.8 - 1)' = na.omit(pa_pAUC_sensitivity) %>% mean(),
                'mean PA pAUC (specificity 0.8 -1)' = na.omit(pa_pAUC_specificity) %>% mean(),
                'mean PA prediction accuracy' = na.omit(pa_prediction_accuracy) %>% mean(),
                'mean PA correlation' = na.omit(pa_correlation) %>% mean(),
                'mean PA kappa' = na.omit(pa_kappa) %>% mean()) %>%
    arrange(-`mean PA AUC`)

  Table1[2:ncol(Table1)] <- Table1[2:ncol(Table1)] %>% round(digits = 3)
  
  Table1 %>%
  write.csv(file = "tables/Table1.csv",
            row.names = FALSE)
  
  
##################################  

  # Table 1 for small sample sizes only
  Table1_sss <-
  full_output %>%
    filter(n_presence <= 20) %>%
    group_by(model) %>%
    summarise('median PA AUC' = na.omit(pa_AUC) %>% median(),
              'mean PA AUC' = na.omit(pa_AUC) %>% mean(),
              'mean PA sensitivity' = na.omit(pa_sensitivity) %>% mean(),
              'mean PA specificity' = na.omit(pa_specificity) %>% mean(),
              'mean PA pAUC (sensitivity 0.8 - 1)' = na.omit(pa_pAUC_sensitivity) %>% mean(),
              'mean PA pAUC (specificity 0.8 -1)' = na.omit(pa_pAUC_specificity) %>% mean(),
              'mean PA prediction accuracy' = na.omit(pa_prediction_accuracy) %>% mean(),
              'mean PA correlation' = na.omit(pa_correlation) %>% mean(),
              'mean PA kappa' = na.omit(pa_kappa) %>% mean()) %>%
    arrange(-`mean PA AUC`)
  
  Table1_sss[2:ncol(Table1_sss)] <- Table1_sss[2:ncol(Table1_sss)] %>% round(digits = 3)
  
  Table1_sss %>%
    write.csv(file = "tables/Table1_sss.csv",
              row.names = FALSE)
  
######################################################################
  
  # Table 2: performance ranked with CV data
  
  Table2 <-
    
    fold_output %>%
    group_by(model) %>%
    summarise('median testing AUC' = na.omit(testing_AUC) %>% median(),
              'mean testing AUC' = na.omit(testing_AUC) %>% mean(),
              'mean testing sensitivity' = na.omit(testing_sensitivity) %>% mean(),
              'mean testing specificity' = na.omit(testing_specificity) %>% mean(),
              'mean testing pAUC (sensitivity 0.8 - 1)' = na.omit(testing_pAUC_sensitivity) %>% mean(),
              'mean testing pAUC (specificity 0.8 -1)' = na.omit(testing_pAUC_specificity) %>% mean(),
              'mean testing prediction accuracy' = na.omit(testing_prediction_accuracy) %>% mean(),
              'mean testing correlation' = na.omit(testing_correlation) %>% mean(),
              'mean testing kappa' = na.omit(testing_kappa) %>% mean()) %>%
    arrange(-`mean testing AUC`)
  
  Table2[2:ncol(Table2)] <- Table2[2:ncol(Table2)] %>% round(digits = 3)
  
  Table2 %>%
    write.csv(file = "tables/Table2.csv",
              row.names = FALSE)
  
  
  
######################################################################
  
  
  
  Table2_sss <-
    fold_output %>%
    filter(n_presence <= 20) %>%
    group_by(model) %>%
    summarise('median testing AUC' = na.omit(testing_AUC) %>% median(),
              'mean testing AUC' = na.omit(testing_AUC) %>% mean(),
              'mean testing sensitivity' = na.omit(testing_sensitivity) %>% mean(),
              'mean testing specificity' = na.omit(testing_specificity) %>% mean(),
              'mean testing pAUC (sensitivity 0.8 - 1)' = na.omit(testing_pAUC_sensitivity) %>% mean(),
              'mean testing pAUC (specificity 0.8 -1)' = na.omit(testing_pAUC_specificity) %>% mean(),
              'mean testing prediction accuracy' = na.omit(testing_prediction_accuracy) %>% mean(),
              'mean testing correlation' = na.omit(testing_correlation) %>% mean(),
              'mean testing kappa' = na.omit(testing_kappa) %>% mean()) %>%
    arrange(-`mean testing AUC`)
  
  Table2_sss[2:ncol(Table2_sss)] <- Table2_sss[2:ncol(Table2_sss)] %>% round(digits = 3)
  
  Table2_sss %>%
    write.csv(file = "tables/Table2_sss.csv",
              row.names = FALSE)
  

######################################################################
  
  # Evaluating on full models
  
  Table3 <-
    
    full_output %>%
    group_by(model) %>%
    summarise('median training AUC' = na.omit(full_AUC) %>% median(),
              'mean training AUC' = na.omit(full_AUC) %>% mean(),
              'mean training sensitivity' = na.omit(full_sensitivity) %>% mean(),
              'mean training specificity' = na.omit(full_specificity) %>% mean(),
              'mean training pAUC (sensitivity 0.8 - 1)' = na.omit(full_pAUC_specificity) %>% mean(),
              'mean training pAUC (specificity 0.8 -1)' = na.omit(full_pAUC_specificity) %>% mean(),
              'mean training prediction accuracy' = na.omit(full_prediction_accuracy) %>% mean(),
              'mean training correlation' = na.omit(full_correlation) %>% mean(),
              'mean training kappa' = na.omit(full_kappa) %>% mean()) %>%
    arrange(-`mean training AUC`)
  
  Table3[2:ncol(Table3)] <- Table3[2:ncol(Table3)] %>% round(digits = 3)
  
  Table3 %>%
    write.csv(file = "tables/Table3.csv",
              row.names = FALSE)
  
##################################  
  
  # Table 3 for small sample sizes only
  Table3_sss <-
    full_output %>%
    filter(n_presence <= 20) %>%
    group_by(model) %>%
    summarise('median training AUC' = na.omit(full_AUC) %>% median(),
              'mean training AUC' = na.omit(full_AUC) %>% mean(),
              'mean training sensitivity' = na.omit(full_sensitivity) %>% mean(),
              'mean training specificity' = na.omit(full_specificity) %>% mean(),
              'mean training pAUC (sensitivity 0.8 - 1)' = na.omit(full_pAUC_specificity) %>% mean(),
              'mean training pAUC (specificity 0.8 -1)' = na.omit(full_pAUC_specificity) %>% mean(),
              'mean training prediction accuracy' = na.omit(full_prediction_accuracy) %>% mean(),
              'mean training correlation' = na.omit(full_correlation) %>% mean(),
              'mean training kappa' = na.omit(full_kappa) %>% mean()) %>%
    arrange(-`mean training AUC`)
  
  Table3_sss[2:ncol(Table3_sss)] <- Table3_sss[2:ncol(Table3_sss)] %>% round(digits = 3)
  
  Table3_sss %>%
    write.csv(file = "tables/Table3_sss.csv",
              row.names = FALSE)
  