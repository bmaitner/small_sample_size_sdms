#' @param presence dataframe of covariates at presence points
#' @param background Optional. Dataframe of covariates at background points
#' @param method Optional. If supplied, both presence and background density estimation will use this method.
#' @param ... Additional parameters passed to internal functions.
#' @details Current methods include:.
#' @export
#' @return List of class "dr_model" containing model objects and metadata needed for projecting the fitted models.
fit_density_ratio <- function(presence = NULL,
                              background = NULL,
                              method = NULL,
                              ...){
  #Check data and method
  
  #Fit the ratio
  dr <- do.call(what = paste('dr_', method, sep = ""),
                list(presence_data = presence,
                     background_data = background,
                     method = "fit",
                     ...))

  #Prepare output
  model <- list(ratio = dr,
                method = method)
  
  class(model) <- "dr_model"
  return(model)

}#End fx
