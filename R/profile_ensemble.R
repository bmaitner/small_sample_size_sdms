
profile_ensemble <- function(csv_file,
                             ensemble = c("kde/kde","rulsif","maxnet"),
                             env,
                             quantile = 0.05,
                             focal_bbox = NULL){
  
  
  # Load csv
  
      occs <- read.csv(csv_file)
    
  # Transform projection  
      
      occs %>%
        st_as_sf(coords=c("lon","lat")) -> occs

      st_crs(occs) <- st_crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
      
      occs <- st_transform(occs,crs = "WGS84")
      
      
  # Check occs against bbox (if supplied)
      
      if(!is.null(focal_bbox)){
        

        suppressWarnings(occs %>%
          st_intersection(focal_bbox %>% st_as_sfc()) -> occs_pruned)
        
          if(nrow(occs_pruned)==0){
            
            output <- data.frame(csv_file = csv_file,
                                 species = basename(csv_file) %>%
                                   gsub(pattern = ".csv",replacement = ""),
                                 n_presences = nrow(occs),
                                 one_vote = NA,
                                 two_votes = NA,
                                 three_votes = NA,
                                 total_votes = NA,
                                 vote_entropy = NA,
                                 note = "Possible non-native")
            
            return(output)
            
          }
        
      }
      
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
                      width = 100000)
    
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
        
        if(i==1){
          out_rast <- prediction_raster
          }else{
            out_rast <- out_rast+prediction_raster
            }


    }# end i loop
    
    
  # Return vector of votes(counts), entropy of votes, #any, #all,#1,#2,#3
  
  
    vals <- data.frame(values(out_rast)) %>%
      filter(votes != 0)

    val_table <- table(vals) %>%
      as.data.frame()

    #ensure all categories in summary table
      
          
      if(nrow(val_table)==0){
        
        val_table <- data.frame(votes = c(1,2,3),Freq = NA)
        
      }else{
        
        
        val_table %>%
          mutate(votes = as.numeric(votes))%>%
          right_join(data.frame(votes = c(1,2,3)),
                     by = "votes")->val_table
        
        
      }
    

    if(is.null(messages)){messages <- NA}        
    

    output <- data.frame(csv_file = csv_file,
                         species = basename(csv_file) %>%
                           gsub(pattern = ".csv",replacement = ""),
                         n_presences = nrow(occs),
                         one_vote = val_table %>% filter(votes == 1) %>% pull(Freq),
                         two_votes =val_table %>% filter(votes == 2) %>% pull(Freq),
                         three_votes = val_table %>% filter(votes == 3) %>% pull(Freq),
                         total_votes = nrow(vals),
                         pred_entropy = Entropy(vals),
                         vote_entropy = Entropy(val_table$Freq),
                         note = messages)
    

  
}
