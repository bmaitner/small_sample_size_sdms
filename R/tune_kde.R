#' @param presences covariates at presence locations
#' @param background covariates at backgruond locations
#' @param selection_criterion to use for ranking
#' @param return_model_info If TRUE, returns a data.frame of model fit info instead of just the best option
tune_kde <- function(presences,
                          background,
                          selection_criterion = "prediction_accuracy",
                          return_model_info = F){
  
  #' @param bwmethod Bandwidth method to use.  One of 'normal-reference' (the default),'cv.ml', or 'cv.ls'
  
  normal_reference_fit <- evaluate_model(presences = presences,
                                         background = background,
                                         method = "kde",
                                         bwmethod = "normal-reference")
  
  cv.ml_fit <- evaluate_model(presences = presences,
                              background = background,
                              method = "kde",
                              bwmethod = "cv.ml")
  
  cv.ls_fit <- evaluate_model(presences = presences,
                              background = background,
                              method = "kde",
                              bwmethod = "cv.ls")
  
  out <- rbind(normal_reference_fit,cv.ml_fit,cv.ls_fit)
  out$bwmethod <- c("normal-reference","cv.ml","cv.ls")
  colnames(out) <- gsub(pattern = "full_",replacement = "",x = colnames(out))
  
  if(!selection_criterion %in% colnames(out)){
    message("selection_criterion doesn't match returned column names, returning NA")
    
    return(data.frame(type = NA))
  }
  
  
  if(!return_model_info){
    
    return(data.frame(type = out$type[which.max(out[,selection_criterion])]))
    
  }
  
  return(out)
  
  
  
}


#evaluate_model()
  #fit_plug_and_play()
    #pnp_kde


cv.ls_fit <- evaluate_model(presences = presences,
                            background = background,
                            method = "kde",
                            bwmethod = "cv.ls")


source("C:/Users/Brian Maitner/Desktop/current_projects/pbsdm/R/pnp_kde.R")



