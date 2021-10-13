#' @name evaluate_model
#' @title Evalute PBSDM range map quality
#' @description This function uses 5-fold, spatially stratified, cross-validation to evaluate distribution model quality.
#' @param presences Presence coordinates in long,lat format.
#' @param background Environmental rasters
#' @param method Optional. If supplied, both presence and background density estimation will use this method.
#' @param presence_method Optional. Method for estimation of presence density.
#' @param background_method Optional. Method for estimation of background density.
#' @param bootstrap Character.  One of "none" (the default, no bootstrapping),
#' "numbag" (presence function is bootstrapped),
#' or "doublebag" (presence and background functions are bootstrapped).
#' @param bootstrap_reps Integer.  Number of bootstrap replicates to use (default is 100)
#' @param quantile Quantile to use for thresholding.  Default is 0.05 (5 pct training presence). Set to 0 for minimum trainin presence (MTP).
#' @param ... Additional parameters passed to internal functions.
#' @note Either `method` or both `presence_method` and `background_method` must be supplied.
#' @details Current plug-and-play methods include: "gaussian", "kde","vine","rangebagging", "lobagoc", and "none".
#' Current density ratio methods include: "ulsif", "rulsif".
#' @importFrom pROC roc auc
#' @export
evaluate_model <- function(presences,
                           background,
                           method = NULL,
                           presence_method = NULL,
                           background_method = NULL,
                           bootstrap = "none",
                           bootstrap_reps = 100,
                           quantile = 0.05,
                           ...){
  
  #Little internal function to handle nulls in method
  robust_in <- function(element,set){
    if(is.null(element)){
      return(FALSE)
    }else{
      if(element %in% set){
        return(TRUE)
      }else{return(FALSE)}
    }
  }#end robust in
  
  
  
  # Check that methods were supplied
  if(is.null(method) & (is.null(presence_method) &
                        is.null(background_method))) {
    stop("Please supply either (1) method, or (2) both presence_method and background_method")
  }
  
  
  # Assign methods if needed
  if(!is.null(method)) {
    
    presence_method <- method
    background_method <- method
    
  }
  

  #Make empty output
  
  out_full <- data.frame(full_AUC = NA,
                         full_pAUC_specificity = NA,
                         full_pAUC_sensitivity = NA,
                         full_correlation = NA,
                         full_sensitivity = NA,
                         full_specificity = NA,
                         full_prediction_accuracy = NA,
                         runtime = NA,
                         full_DOR = NA,
                         full_kappa = NA)
    
    #Fit models
  
    model <- NULL
    
    #If density ratio was supplied
    if(robust_in(element = method,set = c("ulsif","rulsif","maxnet"))){
      
      runtime <- proc.time()
      
      
      tryCatch(expr = model <- fit_density_ratio(presence = presences,
                                                 background = background,
                                                 method = method,
                                                 ...),
               error = function(e){
                 message("problem fitting")
                 return(NULL)
                 }
              )
      
      
      
      runtime <- proc.time() - runtime
      
    }else{
      
      runtime <- proc.time()
      
      tryCatch(expr = model <- fit_plug_and_play(presence = presences,
                                                 background = background,
                                                 method = method,
                                                 presence_method = presence_method,
                                                 background_method = background_method,
                                                 bootstrap = bootstrap,
                                                 bootstrap_reps = bootstrap_reps,
                                                 ...),
               error = function(e){
                 message("problem fitting")
                 return(NULL)
               }
                 
                 
                 )
      
      
      runtime <- proc.time() - runtime
# 
#       model <- fit_plug_and_play(presence = presences,
#                                  background = background,
#                                  method = method,
#                                  presence_method = presence_method,
#                                  background_method = background_method,
#                                  bootstrap = bootstrap,
#                                  bootstrap_reps = bootstrap_reps,bwmethod=bwmethod)
      
    }
    
    
    #If model couldn't be fit, return NAs
    if(is.null(model)){
      
      out_full$runtime = runtime[3]
      return(out_full)
      
    }
    
    
    #Project model to background points
    
    #if(method %in% c("ulsif", "rulsif")){
    if(robust_in(element = method,set = c("ulsif","rulsif","maxnet"))){
      
      predictions <- project_density_ratio(dr_model = model,
                                           data = rbind(presences,background))
      
    }else{
      
      
      predictions <- project_plug_and_play(pnp_model = model,
                                           data = rbind(presences,background))
      
      
    }
    
    #Convert predictions to a raster
    
  
  
    
    suitability_v_occurrence <-
      data.frame(suitability = predictions,
                 occurrence = c(rep(1,nrow(presences)),
                                rep(0,nrow(background))
                                )
                 )
    
    
    #then we feed the suitability and presence/pseudoabscence data into the AUC function to get an AUC
    
    #out$AUC[i] <- get_auc(fold_suitability_v_occurrence = fold_suitability_v_occurrence)$AUC
    
    #Training data
    
    roc_obj <- roc(response =suitability_v_occurrence$occurrence,
                            predictor = suitability_v_occurrence$suitability)
    
    out_full$full_AUC <-roc_obj$auc
    
    
    
    out_full$full_pAUC_specificity <- auc(roc = roc_obj,
                                            partial.auc = c(.8, 1),
                                            partial.auc.correct = TRUE,
                                            partial.auc.focus = "specificity")[[1]]
    
    out_full$full_pAUC_sensitivity <- auc(roc = roc_obj,
                                            partial.auc = c(.8, 1),
                                            partial.auc.correct = TRUE,
                                            partial.auc.focus = "sensitivity")[[1]]
    
    #Thresholded data
    suitability_v_occurrence$thresholded <- 
      threshold_predictions(predictions = suitability_v_occurrence$suitability,
                          presences = suitability_v_occurrence$occurrence)
  
    
    
      
    TP <- length(which(suitability_v_occurrence$occurrence == 1 &
                         suitability_v_occurrence$thresholded == 1))

    FN <- length(which(suitability_v_occurrence$occurrence == 1 &
                         suitability_v_occurrence$thresholded == 0))
    
    TN <- length(which(suitability_v_occurrence$occurrence == 0 &
                         suitability_v_occurrence$thresholded == 0))
    
    FP <- length(which(suitability_v_occurrence$occurrence == 0 &
                         suitability_v_occurrence$thresholded == 1))
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (FP + TN)
    #precision <- TP / (TP + FP)
    DOR <- (TP*TN)/(FP*FN)
    #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
    prediction_accuracy <- (TP+TN)/(TP+TN+FP+FN)
    P_o <- (TP+TN)/(TP+TN+FP+FN)
    Ppres <- ((TP+FP)/(TP+TN+FP+FN))*((TP+FN)/(TP+TN+FP+FN))
    Pabs <- ((FN+TN)/(TP+TN+FP+FN))*((FP+TN)/(TP+TN+FP+FN))
    P_e <- Ppres+Pabs
    kappa <- (P_o - P_e)/(1-P_e)
    
    
    out_full$full_DOR <- DOR
    out_full$full_prediction_accuracy <- prediction_accuracy
    out_full$full_sensitivity <- sensitivity
    out_full$full_specificity <- specificity
    out_full$full_kappa <- kappa
    
    suitability_v_occurrence <- na.omit(suitability_v_occurrence)
    
    out_full$full_correlation <- cor(x = suitability_v_occurrence$suitability,
                                     y = suitability_v_occurrence$occurrence,
                                     method = "pearson")
    out_full$runtime = runtime[3]
    
    
    
    
  return(out_full)
  
  
  
}#end fx
