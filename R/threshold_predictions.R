#code to threshold non-raster data

threshold_predictions <- function(predictions, presences,quantile = 0.05){

  if(length(predictions)!=length(presences)){
    stop("predictions and presences should have equal length")
    }
  
  predictions_at_occurrences <- predictions[which(presences==1)]
  

  threshold <- stats::quantile(x = predictions_at_occurrences,
                               probs = quantile,
                               na.rm = T)  
    
    
  thresholded_predictions <- rep(0,length(predictions))
  thresholded_predictions[which(predictions > threshold)] <- 1
  
  return(thresholded_predictions)  
}


