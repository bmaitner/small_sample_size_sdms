source("R/evaluate_model.R")
source("R/threshold_predictions.R")
library(pbsdm)
library(pROC)

#' @param presences covariates at presence locations
#' @param background covariates at backgruond locations
#' @param selection_criterion to use for ranking
#' @param proportion_vector Vector of proportions between 0 and 1 over which models will be evaluated
#' @param return_model_info If TRUE, returns a data.frame of model fit info instead of just the best option
#' @note Assumes presence-only fitting
#' @note Since rangebagging is implemented as min < x < max, using a fraction of 1 will make all observations "unsuitable"
tune_rangebagging <- function(presences,
                          background,
                          selection_criterion = "prediction_accuracy",
                          proportion_vector = seq(0.1, 0.9, 0.1),
                          return_model_info = F){

  out <- NULL

  for(i in proportion_vector){
    
    fit_i <- evaluate_model(presences = presences,
                            background = background,
                            presence_method = "rangebagging",
                            background_method = "none",
                            p = i)
    
    #Attach p (proportion sampled)
    fit_i$proportion = i
    
    #Update output
    out <- rbind(out,fit_i)
    

  }
  
    
  colnames(out) <- gsub(pattern = "full_",replacement = "",x = colnames(out))
  
  if(!selection_criterion %in% colnames(out)){
    message("selection_criterion doesn't match returned column names, returning NA")
    
    return(data.frame(type = NA))
  }
  
  
  if(!return_model_info){
    
    return(data.frame(type = out$proportion[which.max(out[,selection_criterion])]))
    
  }
  
  return(out)
  
  
  
}


#############################

