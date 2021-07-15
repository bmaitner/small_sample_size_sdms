#' @param method one of either "fit" or "predict"
#' @param object fitted object returned by a pnp_... function. Only needed when method = "predict" 
#' @keywords internal
pnp_none <- function(method, object = NULL, ...){
  
  #Code to check inputs
  
  #Code for fitting
  if(method == "fit"){
  
    model=list(method = "none")
    
    class(model) <- "pnp_estimate"
    return(model)
    
  }
  
  #Code for predicting
  
  if(method == "predict"){
    
    return(invisible(NULL))  
  }
  
}
