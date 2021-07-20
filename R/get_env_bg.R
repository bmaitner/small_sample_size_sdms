#' @param coords Coordinates (long,lat) to extract values for
#' @param env Environmental rasterstack in any projection 
#' @param method Methods for getting bg points. Current option is buffer
#' @param width Numeric.  Width (meters or map units) of buffer
get_env_bg <- function(coords, env, method = "buffer", width) {
  
  #check for bad coords
  
  if(max(coords[,1]) > 180 | min(coords[,1]) < -180){
    message("Problematic coords")
  }
  
  if(max(coords[,2]) > 90 | min(coords[,2]) < -90){
    message("Problematic coords")
  }
  
  coords <- sp::SpatialPoints(coords = coords,
                              proj4string = CRS(projargs = "EPSG:4326"))
  
  coords <- spTransform(x = coords,CRSobj = env@crs)    
  
  #make buffer
  buff <-
  buffer(x = coords,
         width = width)
  
  #remove any partial NAs from buffer(since we don't want to use them)
  buff_rast <- rasterize(y = env,x = buff)
  buffer_cells <- which(getValues(buff_rast)==1)
  
  env <- do.call(rbind,extract(y = buff,x = env))
  
  na_or_not <-
    apply(X = env,
          MARGIN = 1,
          FUN = function(x){
            any(is.na(x))
            
          }
    )
  
  buffer_cells <- buffer_cells[which(!na_or_not)]
  
  env <- env[which(!na_or_not),]
  

  if(dim(env)[1]!=length(buffer_cells)){stop("Something wrong with get_env_bg")}
  
  
  return(test <- list(env = env,
                      bg_cells = buffer_cells))
  
  
}# end fx
