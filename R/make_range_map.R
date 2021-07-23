#inputs: 
  #occurrence records
  #environmental layers
  #settings for modeling
  #setting for thresholding
  #CV options (Random or spatial)

#outputs
  #binary map
  #model (for future projections, optional)
  #some metric of model fit


occurrences <- occs[c("longitude","latitude")]
#' @param occurrences Presence coordinates in long,lat format.
#' @param env Environmental rasters
#' @param method Optional. If supplied, both presence and background density estimation will use this method.
#' @param presence_method Optional. Method for estimation of presence density.
#' @param background_method Optional. Method for estimation of background density.
#' @param bootstrap Character.  One of "none" (the default, no bootstrapping),
#' "numbag" (presence function is bootstrapped),
#' or "doublebag" (presence and background functions are bootstrapped).
#' @param bootstrap_reps Integer.  Number of bootstrap replicates to use (default is 100)
#' @param quantile Quantile to use for thresholding.  Default is 0.05 (5 pct training presence). Set to 0 for minimum trainin presence (MTP).
#' @param ... Additional parameters passed to internal functions.
#' @note Either `method` or both `presence_method` and `background_method` must be supplied.
#' @details Current plug-and-play methods include: "gaussian", "kde","vine","rangebagging", "lobagoc", and "none".
#' Current density ratio methods include: "ulsif", "rulsif".
make_range_map(occurrences,
               env,
               method = NULL,
               presence_method = NULL,
               background_method = NULL,
               bootstrap = "none",
               bootstrap_reps = 100,
               quantile = 0.05){
  
  
  #Get presence and background data
    presence_data <- get_env_pres(coords = occurrences,
                                  env = env)
    
    bg_data <- get_env_bg(coords = occurrences,
                          env = env,
                          method = "buffer",
                          width = NULL)
  
  #If density ratio was supplied
    if(method %in% c("ulsif", "rulsif")){
      
      model <- fit_density_ratio(presence = presence_data$env,
                        background = bg_data$env,
                        method = method)
  
    }else{
      
      
      model <- fit_plug_and_play(presence = presence_data$env,
                                 background = bg_data$env,
                                 method = method,
                                 presence_method = presence_method,
                                 background_method = background_method,
                                 bootstrap = bootstrap,
                                 bootstrap_reps = bootstrap_reps)

      
    }  
    
    
    #Project model to background points
    
    
    if(method %in% c("ulsif", "rulsif")){
      
      predictions <- project_density_ratio(dr_model = model,
                                     data = bg_data$env)
      
    }else{
      
      
      predictions <- project_plug_and_play(pnp_model = model,
                                           data = bg_data$env)
      
      
    }  
    
    #Convert predictions to a raster
      prediction_raster <- setValues(env[[1]],
                          values = NA)
      prediction_raster[bg_data$bg_cells] <- predictions
      
    #Apply thresholding
      prediction_raster <- sdm_threshold(prediction_raster = prediction_raster,
                                         occurrence_sp = presence_data$occurrence_sp,
                                         quantile = quantile)
    
    
  
  
  
}