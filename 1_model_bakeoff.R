library(disdat)
source("R/evaluate_disdat.R")


#Select pnp modules to consider (as both numerator and denominator)
  pnp_components <- c("rangebagging",
                      "gaussian",
                      "kde")
  
# Make full set of hybrid models to consider
  models_to_evaluate <- NULL
  
  for(i in pnp_components){
  for(j in pnp_components){
    models_to_evaluate <- 
      rbind(models_to_evaluate,
          data.frame(presence_method = i,
                     background_method = j)) 

  }}


# Iterate through models
  
full_model_outputs <- NULL
fold_model_output <- NULL
  
  for(i in 1:nrow(models_to_evaluate)){
    
    model_i <- 
    evaluate_disdat(presence_method = models_to_evaluate$presence_method[i],
                    background_method = models_to_evaluate$background_method[i])
    
  
    full_model_outputs <- rbind(full_model_outputs,
                                data.frame(pres_method = models_to_evaluate$presence_method[i],
                                           bg_method = models_to_evaluate$background_method[i],
                                           model_i$full_model_stats))
    
    fold_model_outputs <- rbind(fold_model_outputs,
                                data.frame(pres_method = models_to_evaluate$presence_method[i],
                                           bg_method = models_to_evaluate$background_method[i],
                                           model_i$fold_model_stats))
    
    
  }
  

#ggplots of model stats (facet grid of numerator and denominator)  
  
  
  
  
  