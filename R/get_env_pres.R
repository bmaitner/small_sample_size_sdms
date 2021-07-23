#' @param coords Coordinates (long,lat) to extract values for
#' @param env Environmental rasterstack in any projection 
get_env_pres <- function(coords, env) {

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

return(test <- list(env = extract(y = coords,x = env),
                    occurrence_sp = coords))

  
}# end fx

