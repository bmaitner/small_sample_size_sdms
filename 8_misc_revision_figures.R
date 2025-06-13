# This file contains additional figures that were added during revision



#########################################

# Example map of sites
library(sf)
library(tidyverse)
library(disdat)

regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
data_i <- disData(region = region_i)

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


#Skipping quantile 0.1 model gaussian/none already done


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
                    select(species,pa_sensitivity,full_specificity,pa_AUC)%>%
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
    
  }#pnp
  
  for(j in 1:length(dr_models_to_evaluate)){
    
    
    if(quantile_q %in% quantile_variation_output$quantile &
       dr_models_to_evaluate[j] %in% quantile_variation_output$model
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
                out_i$full_model_stats %>%
                  select(species,pa_sensitivity,full_specificity,pa_AUC)%>%
                  mutate(model = dr_models_to_evaluate[j],
                         quantile = quantile_q) %>%
                  mutate(pa_AUC = as.numeric(pa_AUC)))
    
    
    # save temporary bits
    
    quantile_variation_output %>%
      saveRDS(file = file.path(quantile_tempfile_full))
    
    
  }#dr
  
  
  # save temporary bits
  
    quantile_variation_output %>%
      saveRDS(file = file.path(quantile_tempfile_full))

} #q loop




# Plot
  # sens vs spec, color by model, facet by quantile
  # alternatively, can calculate rank, then plot median rank +/- max/min

################################################################################

  # Comparison with Valavi et al 2021



