#' @param pnp_model A fitted plug-and-play model produced by `fit_plug_and_play`
#' @param data covariate data
#' @return A vector of relative suitabilities evaluates at the covariates supplied in the data object.
project_plug_and_play <- function(pnp_model, data) {

  #Check that pnp_model is the correct class
  if(class(pnp_model) != "pnp_model") {
    stop("Invalid pnp_model supplied.")
  }

  f1_est <- do.call(what = paste('pnp_',pnp_model$f1_method,sep = ""),
                    list(data = data,
                         method = "predict",
                         object = pnp_model$f1))
  
  f0_est <- do.call(what = paste('pnp_',pnp_model$f0_method,sep = ""),
                    list(data = data,
                         method = "predict",
                         object = pnp_model$f0))

  return(S=exp(f1_est - f0_est))
}
