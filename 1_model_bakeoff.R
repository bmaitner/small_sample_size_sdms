library(disdat)
library(np)
library(AUC)
library(rvinecopulib)
library(pROC)
source("R/get_env_pres.R")
source("R/get_env_bg.R")
source("R/fit_plug_and_play.R")
source("R/pnp_gaussian.R")
source("R/pnp_rangebagging.R")
source("R/pnp_none.R")
source("R/pnp_kde.R")
source("R/project_plug_and_play.R")
source("R/fit_density_ratio.R")
source("R/project_density_ratio.R")
source("R/dr_ulsif.R")
source("R/evaluate_range_map.R")
source("R/get_auc.R")
source("R/pnp_gaussian.R")
source("R/pnp_lobagoc.R")
source("R/pnp_vine.R")
source("R/evaluate_disdat.R")
source("R/stratify_spatial.R")


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
fold_model_outputs <- NULL
  
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
  
# save outputs as an RDS object (since it takes so long to re-run)
saveRDS(object = full_model_outputs,
        file = "outputs/bake_off_pnp_full_model_outputs.RDS")

saveRDS(object = fold_model_outputs,
        file = "outputs/bake_off_pnp_fold_model_outputs.RDS")

# ggplots of model stats (facet grid of numerator and denominator)  
  
  
  
  
  