library(tidyverse)
library(lubridate)
quick_and_dirty_ensemble_map <- function(env,
                                         species,
                                         ensemble,
                                         buffer_width = 100000,
                                         quantile = 0.05,
                                         oldest_year = 1970){
  
  # Load csv
  
  occs <- BIEN::BIEN_occurrence_species(species = species)
  
  # Filter
  
  occs <-
  occs %>%
    mutate(year_collect = year(date_collected)) %>%
    filter(date_collected >= oldest_year) %>%
    select(scrubbed_species_binomial,latitude,longitude) %>%
    distinct()
  
  
  
  # Transform projection  
  
  occs %>%
    st_as_sf(coords=c("longitude","latitude")) -> occs
  
  st_crs(occs) <- st_crs("WGS84")
  
  # Convert format for simplicity
  
  coords <- st_coordinates(occs)
  
  occs %>%
    st_drop_geometry()%>%
    mutate(lon = coords[,"X"],
           lat = coords[,"Y"]) -> occs
  
  
  # Get background
  
  bg <- S4DM::get_env_bg(coords = occs[c("lon","lat")],
                          env = env,
                          standardize = FALSE,
                          width = buffer_width)
  
  # Get presence  
  
  pres <- S4DM::get_env_pres(coords = occs[c("lon","lat")],
                              env = env)
  
  # Make message vector
  
  messages <- NULL
  
  # Model vector
  
  for( i in 1:length(ensemble)){
    
    model_i <- ensemble[i]
    
    model_i <- strsplit(x = model_i,split = "/") %>%
      unlist()
    
    # fit model
    if(length(model_i) == 2){
      
      fitted_i <- tryCatch(S4DM::fit_plug_and_play(presence = pres$env,
                                                    background = bg$env,
                                                    presence_method =  model_i[1],
                                                    background_method = model_i[2]),
                           error = function(e){e})
      
    }else{
      
      fitted_i <- tryCatch(S4DM::fit_density_ratio(presence = pres$env,
                                                    background = bg$env,
                                                    method = model_i[1]),
                           error = function(e){e})
      
    }
    
    # if error, move to next
    
    if(inherits(fitted_i,"error")){
      
      messages <- paste(messages, "model",i, "failed to fit;")
      
      prediction_raster <- env[[1]]
      values(prediction_raster) <- 0
      
      names(prediction_raster) <- "votes"
      varnames(prediction_raster) <- "votes"
      
      if(i==1){
        out_rast <- prediction_raster
      }else{
        out_rast <- out_rast+prediction_raster
      }
      
      next
      
    }# end model fitting error handling
    
    # project model
    
    if(length(model_i) == 2){
      
      
      projected_i <- S4DM::project_plug_and_play(pnp_model = fitted_i,
                                                  data = bg$env)
      
      
    }else{
      
      projected_i <- S4DM::project_density_ratio(dr_model = fitted_i,
                                                  data = bg$env)
      
    }
    
    # Make raster
    
    prediction_raster <- env[[1]]
    values(prediction_raster) <- NA
    
    prediction_raster[bg$bg_cells] <- projected_i
    
    names(prediction_raster) <- "votes"
    varnames(prediction_raster) <- "votes"
    
    # Threshold
    
    prediction_raster <- sdm_threshold(prediction_raster = prediction_raster,
                                       occurrence_sf = pres$occurrence_sf,
                                       quantile = quantile)
    
    
    prediction_raster[is.na(prediction_raster)] <- 0
    
    if(i == 1){
      out_rast <- prediction_raster
    }else{
      out_rast <- out_rast+prediction_raster
    }
    
    
  }# end i loop
  
  
  # make map
  
  # Get reasonable box for plotting
  
    world <- st_as_sf(maps::map("world", fill=TRUE, plot =FALSE)) %>%
      st_transform(crs = "WGS84")
  
    fl_bbox_sf <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) %>%
      dplyr::filter(ID == "florida") %>%
      st_transform(crs = "WGS84")%>%
      st_bbox(fl) %>%
      st_as_sfc()

    raster_sf <-out_rast %>%
      terra::subst(from = 0, to = NA)%>%
    terra::as.polygons() %>%
      st_as_sf()%>%
      st_transform(crs = "WGS84")
    
    # convert votes to ordered factor
    
    raster_sf %>%
      mutate(votes = factor(x=votes,
                            levels = c(1,2,3)))->raster_sf

    # total_bbox <- fl_bbox_sf %>%
    #   st_union(raster_sf) %>%
    #   st_buffer(dist = 100000)%>%
    #   st_bbox()
    
    total_bbox <- raster_sf %>%
      #st_union(fl_bbox) %>%
      st_buffer(dist = 100000)%>%
      st_bbox()
    
    
library(tidyterra)
    
    plot_out <-
      ggplot(data = world)+
      geom_sf(data = world)+
      geom_sf(data = raster_sf,mapping = aes(fill=votes))+
      # geom_spatraster(data = out_rast%>%
      #                   terra::subst(from = 0, to = NA),
      #                 na.rm = TRUE,
      #                 maxcell = 1e10)+
      # scale_fill_gradient(na.value=NA)+
      coord_sf(ylim=c(total_bbox[[2]], total_bbox[[4]]),
                      xlim=c(total_bbox[[1]],total_bbox[[3]]))+
      theme_bw()


  return(plot_out)
  
}
