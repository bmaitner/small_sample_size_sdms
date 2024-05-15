# Example ensembles


# Load libraries

library(BIEN)
library(geodata)
library(ggplot2)
library(pbsdm)
library(sf)
library(terra)
library(tidyterra)

# Make a temporary directory to store climate data

  temp <- tempdir()

# Get environmental data

env <- worldclim_global(var = "bio",
                        res = 5,
                        path = temp)

# # Function for occurrence counting
#   
#   get_n_occs <- function(file_list){
#     
#     for(i in 1:length(file_list)){
#     
#         print(i)
#       
#         file_i <- read.csv(file_list[i])
#       
#       # assign species if it isn't part of the metadata
#       
#         if(is.null(file_i$sp)){file_i$sp <- gsub(pattern = ".csv",
#                                                  replacement = "",
#                                                  x = basename(file_list[i]))}
#   
#       # Get the count and other metadata  
#   
#         out_i <- data.frame(file = file_list[i],
#                             species=file_i$sp[1],
#                             n_occs = nrow(file_i))
#       # combine output
#         
#         if(i == 1){
#           out <- out_i
#         }else{
#           out <- rbind(out,out_i)
#         }
#       
#     } # i loop end
#     
#     return(out)
#     
#   } #end n_occs fx
# 
#   # Get info on potential species. Using a dump of small-ranged species from BIEN
#   
#   point_counts <- list.files("data/manual_downloads/BIEN_occs/Users/ctg/Documents/SDMs/BIEN_1123/_inputs/SpeciesCSVs/Points/",
#                              pattern = ".csv",full.names = TRUE) |> get_n_occs()
#   
#   # Empty folder
#   # point_counts_from_rangebag <- list.files("data/manual_downloads/BIEN_occs/Users/ctg/Documents/SDMs/BIEN_1123/_inputs/SpeciesCSVs/PointsFromRangeBag/",
#   #                            pattern = ".csv",full.names = TRUE) |> get_n_occs()
#   
#   rangebag_counts <- list.files("data/manual_downloads/BIEN_occs/Users/ctg/Documents/SDMs/BIEN_1123/_inputs/SpeciesCSVs/RangeBag/",
#                              pattern = ".csv",full.names = TRUE) |> get_n_occs()
#   
#   rangebag_counts_from_ppm <- list.files("data/manual_downloads/BIEN_occs/Users/ctg/Documents/SDMs/BIEN_1123/_inputs/SpeciesCSVs/RangeBagFromPPM/",
#                                 pattern = ".csv",full.names = TRUE) |> get_n_occs()
#   
#   ppm1530_counts <- list.files("data/manual_downloads/BIEN_occs/Users/ctg/Documents/SDMs/BIEN_1123/_inputs/SpeciesCSVs/PPM15_30/",
#                                          pattern = ".csv",full.names = TRUE) |> get_n_occs()
#   
#   ppm_counts  <- list.files("data/manual_downloads/BIEN_occs/Users/ctg/Documents/SDMs/BIEN_1123/_inputs/SpeciesCSVs/PPM/",
#                             pattern = ".csv",full.names = TRUE) |> get_n_occs()
#   
# 
#   all_counts <- rbind(point_counts,
#         rangebag_counts,
#         rangebag_counts_from_ppm,
#         ppm1530_counts,
#         ppm_counts)
  
  # saveRDS(object = all_counts,file = "data/manual_downloads/BIEN_occs/occ_counts.RDS")

  all_counts <- readRDS("data/manual_downloads/BIEN_occs/occ_counts.RDS")
  
  ssss_counts <- all_counts %>%
    filter(n_occs <= 20)
  
  
  # 221000 species with 20 or fewer observations!!
  
  # grab a subset from Florida
  
  fl_species <- BIEN_list_state(country = "United States",state = "Florida")
  
  ssss_counts %>%
    mutate(species = gsub(pattern = "_",
                          replacement = " ",
                          x = species)) -> ssss_counts
  
  fl_counts <-
    ssss_counts %>%
    filter(species %in% fl_species$scrubbed_species_binomial)

# To make things a bit faster and easier, we'll limit ourselves to the 2 variables (mean temperature and annual precipitation)

  env <- env[[c(1,12)]] %>% scale()
  
# Define a WGS84 bbox for Florida
  
  fl <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) %>%
    dplyr::filter(ID == "florida")
  
  fl %>% st_transform(crs = "WGS84")->fl
  
  fl_bbox <- st_bbox(fl)

source("R/profile_ensemble.R")  

  for(i in 1:nrow(fl_counts)){
    
    print(i)
    
    out_i <- profile_ensemble(csv_file = fl_counts$file[i],
                     ensemble = c("kde/kde","rulsif","maxnet"),
                     env = env,
                     quantile = 0.05,
                     focal_bbox = fl_bbox)
  
    if(i == 1){
      
      ensemble_profiles <- out_i
      
    }else{
      
      ensemble_profiles <- bind_rows(ensemble_profiles,out_i)
      
    }
    
  }
    
# review ensembles
  
ensemble_profiles %>%
  ggplot(mapping = aes(x=n_presences,y=vote_entropy))+
  geom_point()
  
head(ensemble_profiles)

ensemble_profiles %>%
  filter(!is.na(vote_entropy))->test

# Here, we'll use the same data as before for Trillium vaseyi.

#First, we'll check the background data to look for correlated precictors

tv_bg <- get_env_bg(coords = tv[c("longitude","latitude")],
                    env = env,
                    width = 50000,
                    standardize = TRUE) #note that we used a small set of background points to expedite model fitting

# Filter bg of highly correlated stuff

library(corrplot)

tv_bg$env %>%
  cor()


# toss anything greater than 0.7
  tv_bg$env %>%
    as.data.frame() %>%
    select(- wc2.1_10m_bio_3,
           - wc2.1_10m_bio_4,
           - wc2.1_10m_bio_5,
           - wc2.1_10m_bio_6,
           - wc2.1_10m_bio_7,
           - wc2.1_10m_bio_8,
           - wc2.1_10m_bio_9,
           - wc2.1_10m_bio_10,
           - wc2.1_10m_bio_11,
           #- wc2.1_10m_bio_13,
            - wc2.1_10m_bio_14,
            - wc2.1_10m_bio_15,
           - wc2.1_10m_bio_16,
           - wc2.1_10m_bio_17,
           - wc2.1_10m_bio_18,
           - wc2.1_10m_bio_19) %>%
    cor() %>%colnames()->preds_to_keep
   #%>% corrplot()

  env <- env[[preds_to_keep]]
    
# get new bg

  tv_bg <- get_env_bg(coords = tv[c("longitude","latitude")],
                      env = env,
                      width = 50000,
                      standardize = TRUE) #note that we used a small set of background points to expedite model fitting

# Next, we get the presence data:

  tv_presence <- get_env_pres(coords = tv[c("longitude","latitude")],
                              env = env,
                              env_bg = tv_bg)

# make maps
    
  kde_map <- make_range_map(occurrences = tv[c("longitude","latitude")],
                 env = env,
                 method = "kde",
                 background_buffer_width = 50000)  
  
  maxnet_map <- make_range_map(occurrences = tv[c("longitude","latitude")],
                            env = env,
                            method = "maxnet",
                            background_buffer_width = 50000)
  
  rulsif_map <- make_range_map(occurrences = tv[c("longitude","latitude")],
                            env = env,
                            method = "rulsif",
                            background_buffer_width = 50000)  

  tv_stack <- c(kde_map,maxnet_map,rulsif_map)

# combine layers

  tv_sf <-  tv_stack %>%
    terra::app(fun = function(x){sum(x,na.rm = TRUE)})%>%
    subst(from = 0,to=NA)%>%
      terra::as.polygons()%>%
      st_as_sf()%>%
    rename(votes = lyr.1)%>%
    mutate(votes = as.factor(votes))
           
# make point sf


  tv[c("longitude","latitude")]%>%
    st_as_sf(coords = c("longitude","latitude")) -> tv_points
    
  st_crs(tv_points) <- "wgs84"

           
library(tidyterra)    
  library(ggnewscale)

  ggplot(data = tv_sf)+
    geom_sf()+
  geom_spatraster(data = env$wc2.1_10m_bio_1 %>%
                    crop(y = tv_sf),na.rm = TRUE)+
    scale_fill_viridis_c(na.value = NA,name="Mean Ann. Temp.")+
    theme_bw()+
    new_scale_fill()+
    geom_sf(data = tv_sf,
            mapping = aes(fill = votes))+
      scale_fill_manual(breaks = c(1,2,3),
                        values= c("purple",
                                  "magenta",
                                  "skyblue"))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    geom_sf(data = tv_points,size=.4,alpha=0.5)
  
  
kde_map%>%
#rulsif_map%>%
  #maxnet_map%>%
  terra::as.polygons()%>%
  plot()
  
###################################




# Opuntia corallicola 
# Pilosocereus robinii 

library(todoBIEN)
source("C:/Users/Brian Maitner/Desktop/current_projects/cds/misc/password_and_username.R")
cm<-todoBIEN::BIEN_occurrence_species(species = "Consolea macracantha",
                                  user=user,
                                  password=password) %>%
  select(latitude,longitude) %>%
  unique()

pr<-todoBIEN::BIEN_occurrence_species(species = "Pilosocereus robinii",
                                      user=user,
                                      password=password) %>%
  select(latitude,longitude) %>%
  unique()


env <- worldclim_global(var = "bio",
                        res = 10,
                        path = temp)

bien_output <-pr
quicklook <- function(bien_output,env){
  
  kde_map_bien <- make_range_map(occurrences = bien_output[c("longitude","latitude")],
                            env = env,
                            method = "kde",
                            background_buffer_width = 50000)  
  
  maxnet_map_bien <- make_range_map(occurrences = bien_output[c("longitude","latitude")],
                               env = env,
                               method = "maxnet",
                               background_buffer_width = 50000)
  
  rulsif_map_bien <- make_range_map(occurrences = bien_output[c("longitude","latitude")],
                               env = env,
                               method = "rulsif",
                               background_buffer_width = 50000)
  
  
  bien_sf <-  c(kde_map_bien,maxnet_map_bien,rulsif_map_bien) %>%
    terra::app(fun = function(x){sum(x,na.rm = TRUE)})%>%
    subst(from = 0,to=NA)%>%
    terra::as.polygons()%>%
    st_as_sf()%>%
    rename(votes = lyr.1)%>%
    mutate(votes = as.factor(votes))
  
  
  return(bien_sf)
    
  
}

#





test <- todoBIEN::BIEN_occurrence_species(species = "Deeringothamnus pulchellus",
                                          user=user,
                                          password=password) %>%
  select(latitude,longitude) %>%
  unique()


quicklook(bien_output = test,
          env = env)->out

ggplot(data = out)+
  geom_sf()+
  geom_spatraster(data = env$wc2.1_10m_bio_1 %>%
                    crop(y = st_buffer(out,1000000)),na.rm = TRUE)+
  scale_fill_viridis_c(na.value = NA,name="Mean Ann. Temp.")+
  theme_bw()+
  new_scale_fill()+
  geom_sf(data = out,
          mapping = aes(fill = votes))+
  scale_fill_manual(breaks = c(1,2,3),
                    values= c("purple",
                              "magenta",
                              "skyblue"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand=c(0,0))



