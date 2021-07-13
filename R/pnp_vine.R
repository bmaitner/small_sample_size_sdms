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


pdevine <- kdevine::kdevine(x = pres_env)

library(kdevine)

library(rvinecopulib)

pvine <- vine(data = pres_env)
?vine

  all_vine<-
  dvine(x = all_env_vals,
      vine = pvine)

univine<-vine(data = pres_env[,2])  

univine_predict <- dvine(x = pres_env[,2],
                         vine = univine)


test <- pnp_vine(data = pres_env,
                 method = "fit")


test$f
test_out <- pnp_vine(data = all_env_vals,
                     method = "predict",
                     object = test$f)

?seq
?sequence
max(pres_env[,2])



predvine <- env[[2]]
predvine <- setValues(predvine,NA)

predvine[which(!na_or_not)] <- test_out
plot(predvine)
plot(log(abs(predvine)))
pred2[which(getValues(pred2)>10)] <- 0

plot(pred2,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))
pred2[which(getValues(pred2)>0.01)] <- 1
plot(pred2)


