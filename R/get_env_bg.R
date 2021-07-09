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
  
  buff <-
  buffer(x = coords,
         width = width)
  
  
  return(test <- do.call(rbind,extract(y = buff,x = env)))
  
  
}# end fx
