#' @param data dataframe of covariates
#' @param method one of either "fit" or "predict"
#' @param object fitted object returned by a pnp_... function. Only needed when method = "predict" 
#' @keywords internal
#' @importFrom rvinecopulib vine
#' @importFrom stats fitted
pnp_vine <- function(data, method, object = NULL){
  
  #Code to check inputs
  
  #Code for fitting
  if(method == "fit"){
    
    f <- vine(data = data)
    
    model <- list(f = f,
                  method = "vine")
    
    class(model) <- "pnp_estimate"
    return(model)
    
  }
  
  #Code for predicting
  
  if(method == "predict"){
    
    #log convert for consistency with other functions
    prediction <- log(dvine(x = data,
                            vine = object$f))
    
    
    return(prediction)  
  }
  
}

##############################################
