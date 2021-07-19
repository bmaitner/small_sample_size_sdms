#' @param data dataframe of covariates
#' @param method one of either "fit" or "predict"
#' @param object fitted object returned by a pnp_... function. Only needed when method = "predict" 
#' @param v Positive integer. The Number of votes to use (default is 100)
#' @param nu Numeric. Tuning parameter for nu-svm
#' @param sigma NULL or Numeric.  Tuning parameter of rbf kernel, will estimate if left NULL (default).
#' @details For fitting, an object is not required (and will be ignored). For prediction, parameters v,p,and d are not needed and will be ignored.
#' @import kernlab kvsm
#' @keywords internal
pnp_lobagoc <- function(data, method, object = NULL, v = 100, nu = 0.01, sigma = NULL){
  
  #Code to check inputs
  
  #Code for fitting
  if(method == "fit"){
    
    labels <- c(rep(1, dim(data)[1]))
    models <- list()
    for(i in 1:v){
      data_i <- rbind(data[sample.int(n = dim(data)[1],
                                 replace = TRUE),])
      if(!is.null(sigma)){

        models[[i]] <- ksvm(labels~.,
                            data = cbind(data, labels),
                            type='one-svc',
                            nu=nu,
                            kpar=list(sigma = sigma))
      }
      else{

        models[[i]] <- ksvm(labels~.,
                            data = cbind(data, labels),
                            type='one-svc',
                            nu = nu)
      }
    }
    

    model <- list(lobagoc_models = models,
                  method = "lobagoc")
    
    class(model) <- "pnp_estimate"
    return(model)
    
  }
  
  #Code for predicting
  
  if(method == "predict"){
    
    v <- length(model$lobagoc_models)
    prediction <- data[,1]*0 

    for(i in 1:n.votes){
      prediction <- prediction + predict(model$lobagoc_models[[i]], data)
    }
    
    prediction <- (prediction/v)
    
    return(log(prediction))  
  }
  
}


