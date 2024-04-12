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
  
#############################################################################
  
  # stats on model performance
  
  library(betareg)
  library(glmmTMB)
  
  # set up data for marginal
  
  marginal_stats <-NULL
  

#AUC
    betareg(data = combined_output %>%
              mutate(n_presence = log10(n_presence))%>%
              filter(!pres_method %in% c("vine"),
                     !bg_method %in% c("lobagoc")) %>%
              dplyr::select(pa_AUC,pres_method,bg_method,n_presence)%>%
              na.omit() %>%
              unique(),
            formula = pa_AUC ~ pres_method + bg_method + n_presence
            + n_presence*pres_method
            + pres_method*bg_method,
            #link = "logit",
            type="BR") -> AUC_model
    
    #pseudor2 = ~.15
  
      
    summary(AUC_model)

  
# function to wrangle marginal effects
    
    library(margins)
      get_marg_df <- function(model){
        
          data.frame(response = model$formula[[2]] %>% as.character(),
                     model %>% margins() %>% summary())
      }
      

      marginal_stats <-
        bind_rows(marginal_stats,
          data.frame(get_marg_df(AUC_model),
                     model_r2 = AUC_model$pseudo.r.squared))
      
      
            
#spec (need to use a model that includes 0 and 1)

      glm(data = combined_output %>%
            mutate(n_presence = log10(n_presence))%>%
          rowwise()%>%
          mutate(n_pa_pts_total = sum(n_pa_absence+n_pa_presence))%>%
          filter(!pres_method %in% c("vine"),
                 !bg_method %in% c("lobagoc")) %>%
          dplyr::select(pa_sensitivity,pres_method,bg_method,n_presence,n_pa_pts_total,species)%>%
          na.omit() %>%
          unique(),
        formula = pa_sensitivity ~ pres_method + bg_method + n_presence
        + n_presence*pres_method
        + pres_method*bg_method,
        family = "binomial"
        #,weights = n_pa_pts_total
        ) -> sens_model

      
      marginal_stats <-
        bind_rows(marginal_stats,
          data.frame(get_marg_df(sens_model),
                     model_r2 = pR2(sens_model)['McFadden']))
      
    glm(data = combined_output %>%
          mutate(n_presence = log10(n_presence))%>%
          
          rowwise()%>%
          mutate(n_pa_pts_total = sum(n_pa_absence+n_pa_presence))%>%
          filter(!pres_method %in% c("vine"),
                 !bg_method %in% c("lobagoc")) %>%
          dplyr::select(pa_specificity,pres_method,bg_method,n_presence,n_pa_pts_total)%>%
          na.omit() %>%
          unique(),
        formula = pa_specificity ~ pres_method + bg_method + n_presence
        + n_presence*pres_method
        + pres_method*bg_method,
        family = "binomial"
        #,weights = n_pa_pts_total
        ) -> spec_model
    
        marginal_stats <-
          bind_rows(marginal_stats,
                    data.frame(get_marg_df(spec_model),
                               model_r2 = pR2(spec_model)['McFadden']))
        
        
        # Marginal plots to ask how (in general) they're different 

        

        marginal_stats %>%
          mutate(plot_order = case_when(str_detect(factor,"bg") ~ 3,
                                  str_detect(factor,"pres_meth") ~ 2,
                                  str_detect(factor,"n_pres") ~ 1))%>%
          mutate(factor = gsub(pattern = "n_presence",replacement = "N Presence",x=factor))%>%
          mutate(factor = gsub(pattern = "pres_method",replacement = "Pres:",x=factor))%>%
          mutate(factor = gsub(pattern = "bg_method",replacement = "BG:",x=factor))%>%
          mutate(response = case_when(response == "pa_AUC" ~ "AUC",
                                      response == "pa_sensitivity" ~ "Sensitivity",
                                      response == "pa_specificity" ~ "Specificity")) -> marg_for_gg
        
        
            marg_for_gg %>%
            mutate(factor = factor(factor,levels=marg_for_gg %>%
                                     arrange(plot_order,factor) %>%
                                     pull(factor) %>%
                                     unique()
            ))%>%
          ggplot(mapping = aes(y=AME,x=factor))+
          geom_point()+
          #scale_x_discrete(labels=marg_for_gg$factor_name)+
          geom_errorbar(mapping = aes(ymin=lower,ymax=upper))+
          facet_wrap(~response,ncol = 1,scales = "free_y")+
          geom_hline(yintercept = 0,lty=2)+
          theme_bw()+
          geom_vline(
            aes(xintercept =  stage(factor, after_scale = 1.5)),
            colour = "grey"
          )+
          geom_vline(
            aes(xintercept =  stage(factor, after_scale = 3.5)),
            colour = "grey"
          ) -> marginals_plot
        
            ggsave(plot = marginals_plot, filename = "figures/marginals_plot.svg",
                   width = 10,height = 5,units = "in",dpi = 600)
            
            ggsave(plot = marginals_plot, filename = "figures/marginals_plot.jpg",
                   width = 10,height = 5,units = "in",dpi = 600)
            

        
        # Anova of homotypic vs heterotypic models in general? Or t-test?
        # ANOVA of model (pres+bg as a single value) to test for difference in AUC, etc
          # may be best to import to GGPLOT to make prettier
