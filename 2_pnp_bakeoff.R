# 2a pnp bakeoff analyses

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
    filter(model != "CVmaxnet" ) %>%
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
    filter(model != "CVmaxnet" ) %>%
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
    filter(model != "CVmaxnet" ) %>%
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
    filter(model != "CVmaxnet" ) %>%
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
    filter(model != "CVmaxnet" ) %>%
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
  
    # function to wrangle marginal effects
    
      library(margins)
  
      get_marg_df <- function(model){
        
        data.frame(response = model$formula[[2]] %>% as.character(),
                   model %>% margins() %>% summary())
      }
  
    # function to wrangle interaction effects
      
      source("R/get_int_effects.R")

  
  # plot to identify full-rank
  
  full_output %>%
    filter(model != "CVmaxnet" ) %>%
    ggplot(mapping = aes(x = n_presence,
                         y = pa_AUC))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)
  
  
  library(betareg)
  library(glmmTMB)
  library(pscl)
  
  # set up data for marginal
  
  marginal_stats <-NULL
  int_stats <- NULL  

#AUC
    betareg(data = full_output %>%
              filter(!is.na(pres_method))%>%
              mutate(n_presence = log10(n_presence))%>%
              filter(!bg_method %in% c("lobagoc")) %>%
              dplyr::select(pa_AUC,pres_method,bg_method,n_presence)%>%
              na.omit() %>%
              unique(),
            formula = pa_AUC ~ pres_method + bg_method + n_presence
            + n_presence*pres_method
            + pres_method*bg_method,
            #link = "logit",
            type="BR") -> AUC_model
  
  
    glmmTMB(data = full_output %>%
              filter(!is.na(pres_method))%>%
              mutate(n_presence = log10(n_presence))%>%
              filter(!bg_method %in% c("lobagoc")) %>%
              dplyr::select(pa_AUC,pres_method,bg_method,n_presence)%>%
              na.omit() %>%
              unique(),
            formula = pa_AUC ~ pres_method + bg_method + n_presence
            + n_presence*pres_method
            + pres_method*bg_method,
            family = beta_family()) -> AUC_model_alt

    #pseudor2 = ~.11

    summary(AUC_model)

      marginal_stats <-
        bind_rows(marginal_stats,
          data.frame(get_marg_df(AUC_model),
                     model_r2 = AUC_model$pseudo.r.squared))
      
      int_stats <-
        bind_rows(int_stats,
                  data.frame(get_int_effects(model = AUC_model_alt))
                  )

#spec (need to use a model that includes 0 and 1)

      glm(data = full_output %>%
            filter(!is.na(pres_method))%>%
            mutate(n_presence = log10(n_presence)) %>%
          rowwise()%>%
          mutate(n_pa_pts_total = sum(n_pa_absence+n_pa_presence))%>%
          # filter(!pres_method %in% c("vine", "lobagoc"),
          #        !bg_method %in% c("lobagoc", "vine")) %>%
           dplyr::select(pa_sensitivity,pres_method,bg_method,n_presence,n_pa_pts_total,species)%>%
          na.omit() %>%
          unique(),
        formula = pa_sensitivity ~ pres_method + bg_method + n_presence
        + n_presence*pres_method
        + pres_method*bg_method,
        family = "binomial"
        #,weights = n_pa_pts_total
        ) -> sens_model
      
      summary(sens_model)
      
      marginal_stats <-
        bind_rows(marginal_stats,
          data.frame(get_marg_df(sens_model),
                     model_r2 = pR2(sens_model)['McFadden']))
      
      int_stats <-
        bind_rows(int_stats,
                  data.frame(get_int_effects(model = sens_model))
                  )
      
      
    glm(data = full_output %>%
          mutate(n_presence = log10(n_presence))%>%
          rowwise()%>%
          mutate(n_pa_pts_total = sum(n_pa_absence+n_pa_presence))%>%
          # filter(!pres_method %in% c("vine","lobagoc"),
          #        !bg_method %in% c("lobagoc","vine")) %>%
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

        int_stats <-
          bind_rows(int_stats,
                    data.frame(get_int_effects(model = spec_model))
          )
        
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
            mutate(factor = factor(factor,
                                   levels=marg_for_gg %>%
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
            aes(xintercept =  stage(factor, after_scale = 5.5)),
            colour = "grey"
          ) -> marginals_plot
        
            ggsave(plot = marginals_plot, filename = "figures/marginals_plot.svg",
                   width = 10,height = 5,units = "in",dpi = 600)
            
            ggsave(plot = marginals_plot, filename = "figures/marginals_plot.jpg",
                   width = 10,height = 5,units = "in",dpi = 600)
            
            
      # Interaction plots
            
      
      int_stats %>%
        mutate(model = paste(pres_method,"/",bg_method,sep = ""))%>%
        mutate(response = case_when(response == "pa_AUC" ~ "AUC",
                                    response == "pa_sensitivity" ~ "Sensitivity",
                                    response == "pa_specificity" ~ "Specificity"))%>%
        ggplot(mapping = aes(x=model,
                             y=fit))+
        geom_point()+
        geom_errorbar(mapping = aes(ymin = fit - se.fit,
                                    ymax = fit + se.fit))+
        facet_wrap(~response,
                   ncol = 1,
                   scales = "free_y")+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 6.5)),
          colour = "grey"
        )+
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 12.5)),
          colour = "grey"
        ) +
        
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 18.5)),
          colour = "grey"
        ) +
        geom_vline(
          aes(xintercept =  stage(model, after_scale = 24.5)),
          colour = "grey"
        ) +
        
        
        ylab(NULL) +
        xlab("Model") -> interactions_plot

      
      ggsave(plot = interactions_plot, filename = "figures/interactions_plot.svg",
             width = 10,height = 6,units = "in",dpi = 600)
      
      ggsave(plot = interactions_plot, filename = "figures/interactions_plot.jpg",
             width = 10,height = 6,units = "in",dpi = 600)
      
            
##################
      
      full_output %>%
        filter(n_presence < 20)%>%
        filter(!is.na(pres_method) & !is.na(bg_method))%>%
        mutate(type = case_when(pres_method == bg_method ~ "homotypic",
                                pres_method != bg_method ~ "heterotypic"))%>%
          mutate(type = case_when(bg_method == "none" ~ "presence only",
                                  bg_method != "none" ~ type))%>%
        group_by(model)%>%
          ggplot(mapping = aes(x=type,y=pa_AUC,fill=bg_method))+
        geom_boxplot()+
        #geom_violin()+
        geom_vline(
          aes(xintercept =  stage(type, after_scale = 1.5)),
          colour = "grey"
        )+
        geom_vline(
          aes(xintercept =  stage(type, after_scale = 2.5)),
          colour = "grey"
        )+
        facet_wrap(~pres_method,scales = "free_y",nrow = 1)
        

##################################################################
      
        
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
        write.csv(file = "tables/small_sample_size_comparison_to_maxnet.csv",
                  row.names = FALSE)
      
##################
  
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
      
                    