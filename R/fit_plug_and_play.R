#' @param presence dataframe of covariates at presence points
#' @param background dataframe of covariates at background points
#' @param method Optional. If supplied, both presence and background density estimation will use this method.
#' @param presence_method Optional. Method for estimation of presence density.
#' @param background_method Optional. Method for estimation of background density.
#' @note Either `method` or both `presence_method` and `background_method` must be supplied.
#' @details Current methods include: "gaussian"
#' @export
fit_plug_and_play <- function(presence,
                              background,
                              method = NULL,
                              presence_method = NULL,
                              background_method = NULL){
  
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
  
  
  # Check that methods are available
  
  #for now do this manually, but once function skeleton is working do this by looking up available internals
  current_modules <- c("gaussian")
  
  if(!presence_method %in% current_modules) {
    stop(paste("Presence method not implemented. Please select one of: ",
               current_modules,".", sep = ""))
  }
  
  if(!background_method %in% current_modules) {
    stop(paste("Background method not implemented. Please select one of: ",
               current_modules,".", sep = ""))
  }
  
  # Fit model components
  f1 <- do.call(what = paste('pnp_', method, sep = ""),
                list(data = presence, method = "fit"))
  
  f0 <- do.call(what = paste('pnp_', method, sep = ""),
                list(data = background, method = "fit"))
  
  model <- list(f1 = f1,
                f0 = f0,
                f1_method = presence_method,
                f0_method = background_method)
  
  class(model) <- "pnp_model"
  return(model)
  
  
}#end fx  
