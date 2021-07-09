#' @param data dataframe of covariates
#' @param method one of either "fit" or "predict"
#' @param object fitted object returned by a pnp_... function. Only needed when method = "predict" 
#' @keywords internal
pnp_gaussian <- function(data, method, object = NULL){
  
  #Code to check inputs
  
  #Code for fitting
  if(method == "fit"){
    
    mean <- colMeans(data) # estimated mean of presence points
    sigma <- stats::cov(data) # estimated covariance of presence points
    
    model=list(mean = mean,
               sigma = sigma,
               method = "gaussian")
    
    class(model) <- "pnp_estimate"
    return(model)
  
  }

  #Code for predicting
  
  if(method == "predict"){

    prediction <- mvtnorm::dmvnorm(data,
                                   mean = object$mean,
                                   sigma = object$sigma,
                                   log = TRUE)
    
    
    return(prediction)  
  }

}
