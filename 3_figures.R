library(tidyverse)

tempfile_full <- "outputs/temp_bakeoff_output_full.rds"
tempfile_fold <- "outputs/temp_bakeoff_output_fold.rds"

tempfile_full_dr <- "outputs/temp_bakeoff_output_full_dr.rds"
tempfile_fold_dr <- "outputs/temp_bakeoff_output_fold_dr.rds"

full_output <- readRDS(tempfile_full)
full_output_dr <- readRDS(tempfile_full_dr)

full_output %>%
  mutate(pa_AUC = as.numeric(pa_AUC),
         full_AUC = as.numeric(full_AUC)) -> full_output

full_output_dr %>%
  mutate(pa_AUC = as.numeric(pa_AUC),
         full_AUC = as.numeric(full_AUC)) -> full_output_dr


# full_output_dr %>%
#   rename(pres_method = ratio_method) %>%
#   mutate(bg_method = "none") %>%
#   bind_rows(full_output) -> full_output

full_output_dr %>%
  bind_rows(full_output)->full_output

full_output %>%
  mutate(model = case_when(!is.na(ratio_method) ~ ratio_method,
         is.na(ratio_method) ~ paste(pres_method, "/", bg_method))) -> full_output

#plot entropy vs model type?

library(ggplot2)
library(tidyverse)

# Examine mean entropy (used for setting factor order)

  full_output %>%
  group_by(model)%>%
    filter(!is.nan(entropy))%>%
    summarise(mean_ent = mean(entropy,na.rm = TRUE))%>%
    arrange(mean_ent)

# Convert model type to ordered factor (factor order based on mean entropy, low to high)
  
  full_output <-
  full_output %>%
    mutate(pres_method = factor(pres_method, levels = c("kde","gaussian","rangebagging","lobagoc")),
           bg_method = factor(bg_method, levels = c("kde","gaussian","rangebagging","lobagoc","none")))

# Entropy plot

ent <-  full_output %>%
    filter(!is.na(bg_method)) %>%
    ggplot(mapping = aes(x=log10(n_presence),y=entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")

# Sensitivity (correct presences)
  
sens<-  full_output %>%
  filter(!is.na(bg_method)) %>%
  
    ggplot(mapping = aes(x=log10(n_presence),y=pa_sensitivity,color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("Sensitivity (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()

# Specificity (correct absences)
  
spec <-  full_output %>%
  filter(!is.na(bg_method)) %>%
  
    ggplot(mapping = aes(x=log10(n_presence),y=pa_specificity, color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("Specificity (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()
  


# AUC
  
auc <-  full_output %>%
  filter(!is.na(bg_method)) %>%
  
    ggplot(mapping = aes(x=log10(n_presence),y=pa_AUC,
                         color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    geom_hline(yintercept = .5,lty=2)+
    xlab("log10(presences)")+
    ylab("AUC (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()
  

  
  # Pred accuracy
  
acc <-  full_output %>%
  filter(!is.na(bg_method)) %>%
  
    ggplot(mapping = aes(x=log10(n_presence),y=pa_prediction_accuracy,
                         color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("Prediction Accuracy (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()
  
  
# DOR
  
dor <-  full_output %>%
  filter(!is.na(bg_method)) %>%
  
    ggplot(mapping = aes(x=log10(n_presence),
                         y=log10(pa_DOR),
                         color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("Diagnostic Odds Ratio (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()  
    
  
# runtime
  
run <-  full_output %>%
  filter(!is.na(bg_method)) %>%
  
    ggplot(mapping = aes(x=log10(n_presence),
                         y=log10(runtime),
                         color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("log10(Runtime s)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()    

library(ggpubr)

# combined figure

  ggarrange(sens,spec,auc,run,
            common.legend = TRUE,
            legend = "bottom")

  
####################################
  
  
  # Figure of the different pnp presence/background methods
  
  # Entropy plot
  
  ent_pnp <-  full_output %>%
    filter(!is.na(bg_method)) %>%
    
    filter(!pres_method %in% c("maxnet","ulsif","rulsif")) %>%
    ggplot(mapping = aes(x=log10(n_presence),y=entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")
  
  # Sensitivity (correct presences)
  
  sens_pnp<-  full_output %>%
    filter(!is.na(bg_method)) %>%
    
    filter(!pres_method %in% c("maxnet","ulsif","rulsif")) %>%
    ggplot(mapping = aes(x=log10(n_presence),y=pa_sensitivity,color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("Sensitivity (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()
  
  # Specificity (correct absences)
  
  spec_pnp <-  full_output %>%
    filter(!is.na(bg_method)) %>%
    
    filter(!pres_method %in% c("maxnet","ulsif","rulsif")) %>%
    ggplot(mapping = aes(x=log10(n_presence),y=pa_specificity, color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("Specificity (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()
  
  
  
  # AUC
  
  auc_pnp <-  full_output %>%
    filter(!is.na(bg_method)) %>%
    
    filter(!pres_method %in% c("maxnet","ulsif","rulsif")) %>%
    ggplot(mapping = aes(x=log10(n_presence),y=pa_AUC,
                         color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    geom_hline(yintercept = .5,lty=2)+
    xlab("log10(presences)")+
    ylab("AUC (P/A)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()
  
  # runtime
  
  run_pnp <-  full_output %>%
    filter(!is.na(bg_method)) %>%
    
    filter(!pres_method %in% c("maxnet","ulsif","rulsif")) %>%
    ggplot(mapping = aes(x=log10(n_presence),
                         y=log10(runtime),
                         color= entropy))+
    geom_point()+
    facet_grid(pres_method ~ bg_method)+
    scale_y_continuous(sec.axis =
                         sec_axis(~ . , name = "Presence Method",
                                  labels = NULL, breaks = NULL))+
    scale_x_continuous(sec.axis =
                         sec_axis(~ . , name = "Background Method",
                                  labels = NULL, breaks = NULL))+
    xlab("log10(presences)")+
    ylab("log10(Runtime s)")+
    scale_color_gradient(low = "#2efcff",high="magenta")+
    theme_bw()    
  
  
  #rainbow plot
    ggarrange(sens_pnp,spec_pnp,auc_pnp,run_pnp,
              common.legend = TRUE,
              legend = "bottom")

    
#############
    
    #pnp types, models by color, using splines
    

    # Sensitivity (correct presences)
    
    sens_pnp_line <-  
      full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),y=pa_sensitivity,color= pres_method))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=pres_method),se = FALSE)+
      facet_grid(~ bg_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Background Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Sensitivity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()
    

    spec_pnp_line <-  
    full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),y=pa_specificity,color= pres_method))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=pres_method),se = FALSE)+
      facet_grid(~ bg_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Background Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Specificity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()
    
    auc_pnp_line <-  
    full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),y=pa_AUC,color= pres_method))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=pres_method),se = FALSE)+
      facet_grid(~ bg_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Background Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("AUC (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()

    ent_pnp_line <-  
    full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),y=entropy,color= pres_method))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=pres_method),se = FALSE)+
      facet_grid(~ bg_method) +
      scale_y_continuous(limits = c(0,12),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Background Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Entropy")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()
    
    ggarrange(sens_pnp_line,spec_pnp_line,auc_pnp_line,ent_pnp_line,
              common.legend = TRUE,
              legend = "bottom")
############
    
  # PNP, model-class focused
    
    
    # Sensitivity (correct presences)
    
    sens_pnp_line_alt <-  
      full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=bg_method,
                           y= pa_sensitivity))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=bg_method),se = FALSE)+
      facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Sensitivity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    
    spec_pnp_line_alt <-  
      full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=bg_method,
                           y= pa_specificity))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=bg_method),se = FALSE)+
      facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Specificity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    

    auc_pnp_line_alt <-  
      full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=bg_method,
                           y= pa_AUC))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=bg_method),se = FALSE)+
      facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("AUC (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    
    
    ent_pnp_line_alt <-  
      full_output %>%
      filter(!is.na(bg_method)) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=bg_method,
                           y= entropy))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=bg_method),se = FALSE)+
      facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,10),expand = c(0,.5))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Entropy")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    
    
    ggarrange(sens_pnp_line_alt,
              spec_pnp_line_alt,
              auc_pnp_line_alt,
              ent_pnp_line_alt,
              common.legend = TRUE,
              legend = "bottom")
    
  
### Dr models vs other models    
    
    sens_dr_line_alt <-  
      full_output %>%
      filter(model %in% c("maxnet","ulsif","rulsif",
                          "rangebagging / none",
                          "gaussian / gaussian",
                          "kde / kde")) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=model,
                           y= pa_sensitivity))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE)+
      #facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Sensitivity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")

    spec_dr_line_alt <-  
      full_output %>%
      filter(model %in% c("maxnet","ulsif","rulsif",
                          "rangebagging / none",
                          "gaussian / gaussian",
                          "kde / kde")) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=model,
                           y= pa_specificity))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE)+
      #facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Specificity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    
    auc_dr_line_alt <-  
      full_output %>%
      filter(model %in% c("maxnet","ulsif","rulsif",
                          "rangebagging / none",
                          "gaussian / gaussian",
                          "kde / kde")) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=model,
                           y= pa_AUC))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE)+
      #facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("AUC (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    

    ent_dr_line_alt <-  
      full_output %>%
      filter(model %in% c("maxnet","ulsif","rulsif",
                          "rangebagging / none",
                          "gaussian / gaussian",
                          "kde / kde")) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=model,
                           y= entropy))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE)+
      #facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,11),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Entropy")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    
    runtime_dr_line_alt <-  
      full_output %>%
      filter(model %in% c("maxnet","ulsif","rulsif",
                          "rangebagging / none",
                          "gaussian / gaussian",
                          "kde / kde")) %>%
      ggplot(mapping = aes(x=log10(n_presence),
                           color=model,
                           y= runtime))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE)+
      #facet_grid(~ pres_method) +
      #scale_y_continuous(limits = c(0,11),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_continuous(sec.axis =
                           sec_axis(~ . , name = "Presence Method",
                                    labels = NULL, breaks = NULL))+
      xlab("log10(presences)")+
      ylab("Runtime (s)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Background\nMethod")
    
    
#####################

    # Rarified data
    
    out_rarified_full <- readRDS("outputs/temp_rarified_full.RDS")
    out_rarified_fold <- readRDS("outputs/temp_rarified_fold.RDS")
    
    

    
    
#####################        
    #Useful to plot AUC vs pa AUC? (show overfitting?)
    
    
    sens_rare_line <-  
      out_rarified_full %>%
      drop_na()%>%
      ggplot(mapping = aes(x=n_presence,
                           color=model,
                           y= pa_sensitivity))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE,method = "loess")+
      #facet_grid(~ pres_method) +
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      # scale_y_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      scale_x_sqrt()+
      # scale_x_continuous(sec.axis =
      #                      sec_axis(~ . , name = "Presence Method",
      #                               labels = NULL, breaks = NULL))+
      xlab("Presences")+
      ylab("Sensitivity (P/A)")+
      #scale_color_gradient(low = "#2efcff",high="magenta")+
      theme_bw()+
      labs(color = "Model")
    
    spec_rare_line <-  
      out_rarified_full %>%
      drop_na()%>%
      ggplot(mapping = aes(x=n_presence,
                           color=model,
                           y= pa_specificity))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE,method = "loess")+
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      scale_x_sqrt()+
      xlab("Presences")+
      ylab("Specificity (P/A)")+
      theme_bw()+
      labs(color = "Model")    
    
    auc_rare_line <-  
      out_rarified_full %>%
      drop_na()%>%
      ggplot(mapping = aes(x=n_presence,
                           color=model,
                           y= pa_AUC))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE,method = "loess")+
      scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      scale_x_sqrt()+
      xlab("Presences")+
      ylab("AUC (P/A)")+
      theme_bw()+
      labs(color = "Model")    
    
    runtime_rare_line <-  
      out_rarified_full %>%
      drop_na()%>%
      ggplot(mapping = aes(x=n_presence,
                           color=model,
                           y= runtime))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE,method = "loess")+
      #scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      scale_x_sqrt()+
      xlab("Presences")+
      ylab("Runtime (s)")+
      theme_bw()+
      labs(color = "Model")
    
    entropy_rare_line <-  
      out_rarified_full %>%
      drop_na()%>%
      ggplot(mapping = aes(x=n_presence,
                           color=model,
                           y= entropy))+
      geom_point(alpha=0.1)+
      geom_smooth(aes(color=model),se = FALSE,method = "loess")+
      #scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      scale_x_sqrt()+
      xlab("Presences")+
      ylab("Entropy")+
      theme_bw()+
      labs(color = "Model")
    
    ggarrange(sens_rare_line,
              spec_rare_line,
              auc_rare_line,
              entropy_rare_line,
              common.legend = TRUE,
              legend = "bottom")
    
    
    out_rarified_full %>%
      group_by(model)%>%
      summarise(n=n())
    