#input
  #background suitability vector
  #presence vs abscence vector
    
#output
  #AUC
  #sensitivity and specificity vs a vector of threshold values

#' @param fold_suitability_v_occurrence Dataframe produced in the evaluate_range_maps workflow
#' @param return_AUC_data optionally returns sensitivity, specificity, and the associated thresholds
get_auc <- function(fold_suitability_v_occurrence,return_AUC_data = FALSE){
  
  # log_seq <- seq(from = log(min(fold_suitability_v_occurrence$suitability[which(fold_suitability_v_occurrence$occurrence==1)])),
  #     to = log(max(fold_suitability_v_occurrence$suitability[which(fold_suitability_v_occurrence$occurrence==1)])),
  #     length.out = 100)
  #
  
  log_seq <- seq(from = min(log(fold_suitability_v_occurrence$suitability[which(fold_suitability_v_occurrence$suitability > 0)])),
                 to = log(max(fold_suitability_v_occurrence$suitability)),
                 length.out = 100)

  
  #need to extend log seq below the lowest threshold value
  
  out <- data.frame(threshold = exp(log_seq),
                   sensitivity = NA,
                   specificity = NA)
  
  for( i in 1:length(log_seq)){
    log_suit <- exp(log_seq[i])
    
    #Sensitivity
    
    sensitivity <- 
    
    length(which(fold_suitability_v_occurrence$occurrence == 1 &
                   fold_suitability_v_occurrence$suitability >= log_suit))/ #true positives
      length(which(fold_suitability_v_occurrence$occurrence == 1)) #Number of predicted presences
    
    #Specificity

    specificity <-

      length(which(fold_suitability_v_occurrence$occurrence == 0 &
               fold_suitability_v_occurrence$suitability < log_suit))/ #correct negatives
      length(which(fold_suitability_v_occurrence$occurrence == 0)) #true negatives
      
    
    out$sensitivity[i] <- sensitivity
    out$specificity[i] <- specificity
    
    
    
  }
  
  
  out$one_minus_specificity <- 1 - out$specificity
  
  out <- out[order(out$one_minus_specificity),]
  
  AUC <- 0
  
  for(i in 1:nrow(out)){
    
    if(i == nrow(out)){next}
    
    
    rec_area <- (out$one_minus_specificity[i+1] -
                   out$one_minus_specificity[i]) *
      out$sensitivity[i]
    
    tri_area <- (out$one_minus_specificity[i+1] -
                   out$one_minus_specificity[i]) *
                (out$sensitivity[i+1] -
                   out$sensitivity[i]) *
                0.05
  
    AUC <-AUC + rec_area + tri_area
      
  }
  
  output <- list(AUC = AUC)
  
  if(return_AUC_data == TRUE){
    output <- list(AUC = AUC, sensitivity_and_specificity = out)
  }else{
    output <- list(AUC = AUC)
  }
    
    
  
  return(output)
  
  
  
  
}
