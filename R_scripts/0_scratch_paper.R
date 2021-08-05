

# get a set of species occurrence data

  occurrence_files <-
  list.files(path = "C:/Users/Brian Maitner/Desktop/current_projects/BIENWorkflow/inst/extdata/Demo_SpeciesCSVs/",
             pattern = ".csv",
             full.names = T)

# get a set of climate data
  library(raster)
  env <- stack(x = "C:/Users/Brian Maitner/Desktop/current_projects/BIENWorkflow/inst/extdata/Demo_Env/AllEnv.tif")
  names(env) <- read.csv("C:/Users/Brian Maitner/Desktop/current_projects/BIENWorkflow/inst/extdata/Demo_Env/layerNames.csv",header = F)[,1]
  
# get continent data
  library(sp)
  library(maptools)
  continents <- rgdal::readOGR(dsn = "C:/Users/Brian Maitner/Desktop/current_projects/global_plant_phylo/data/continents/4a7d27e1-84a3-4d6a-b4c2-6b6919f3cf4b202034-1-2zg7ul.ht5ut.shp")
  continents <- spTransform(x = continents,CRSobj = env@crs)
    
# load Drake Lab's plug n play package to check out
  #library(PlugNPlay)
  
  
##################

  #Quick and dirty sampling of presences and absences
  
  library(np)
  library(rvinecopulib)
  source("R/get_env_pres.R")
  source("R/get_env_bg.R")
  source("R/fit_plug_and_play.R")
  source("R/pnp_gaussian.R")
  source("R/pnp_rangebagging.R")
  source("R/pnp_none.R")
  source("R/pnp_kde.R")
  source("R/project_plug_and_play.R")
  source("R/fit_density_ratio.R")
  source("R/project_density_ratio.R")
  source("R/dr_ulsif.R")
  source("R/evaluate_range_map.R")
  source("R/get_auc.R")
  source("R/pnp_gaussian.R")
  source("R/pnp_lobagoc.R")
  source("R/pnp_vine.R")

  #env <- scale(env)
  all_env <- getValues(env)
  all_env_vals <-na.omit(all_env)
  
  na_or_not <-
    apply(X = all_env,
          MARGIN = 1,
          FUN = function(x){
            any(is.na(x))
            
          }
    )
  
  #occs <- read.csv(occurrence_files[3])#some acacia
  occs <- read.csv(occurrence_files[4])#some acaena
  
  pres_env <- get_env_pres(coords = occs[c("longitude","latitude")],
                           env = env)
  
  pres_env$env <- na.omit(pres_env$env)
  
  bg_env <- get_env_bg(coords = occs[c("longitude","latitude")],
                       env = env,
                       constraint_regions = continents) 
  
  
#########  

  #Making example plots of how the methods work 
  gaussian_mod <- pnp_gaussian(data = pres_env$env[,c(2,4),drop=FALSE],
                      method = "fit")

  kde_mod <- pnp_kde(data = pres_env$env[,c(2,4),drop=FALSE],
                     method = "fit")
  
  bag_mod <- pnp_rangebagging(data = pres_env$env[,c(2,4),drop=FALSE],
                     method = "fit")
  
  lobagoc_mod <- pnp_lobagoc(data = pres_env$env[,c(2,4),drop=FALSE],
                              method = "fit")
  
  vine_mod <- pnp_vine(data = pres_env$env[,c(2,4),drop=FALSE],
                              method = "fit")
  sample_data_bivariate <- NULL
  for(i in seq(from=0,to=40,by=1)){
    for(j in seq(from = 0,
                 to = 4000,
                 by = 100)){
      
      sample_data_bivariate <- rbind(sample_data_bivariate,cbind(i,j))
    }
    
  }
  
  
  # 2=0:40
  # 4=400:4000
  
  sample_data <- seq(from = 0, to = 40, by = .1)
  sample_data <- as.matrix(sample_data)
  
  colnames(sample_data) <- colnames(pres_env$env[,2,drop=FALSE])
  colnames(sample_data_bivariate) <- colnames(pres_env$env[,c(2,4)])
  
  gaussian_data <- pnp_gaussian(data = sample_data_bivariate,method = "predict",object = gaussian_mod)
  kde_data <- pnp_kde(data = sample_data_bivariate,method = "predict",object = kde_mod)
  bag_data <- pnp_rangebagging(data = sample_data_bivariate,method = "predict",object = bag_mod)
  lobagoc_data <- pnp_lobagoc(data = sample_data_bivariate,method = "predict",object = lobagoc_mod)
  vine_data <- pnp_vine(data = sample_data_bivariate,
                        method = "predict",
                        object = vine_mod) #fix
  
  combined_data <- cbind(sample_data_bivariate,
                         gaussian_data,
                         kde_data,
                         bag_data,
                         lobagoc_data,
                         vine_data)
  colnames(combined_data)[which(colnames(combined_data)=="")] <- "lobagoc"
  colnames(combined_data) <- gsub(pattern = "_data",
                                  replacement = "",
                                  x = colnames(combined_data))
  
  combined_data <- as.data.frame(combined_data)
  combined_data[3:ncol(combined_data)] <- exp(combined_data[3:ncol(combined_data)])
  
  #combined_data[3:ncol(combined_data)] <- scale(combined_data[3:ncol(combined_data)])
  
  
  library(tidyverse)
  
  combined_data %>%
    pivot_longer(cols = 3:ncol(combined_data)) %>%
    mutate(value = exp(value)) -> combined_data
  
  library(ggplot2)
  
  ggplot(data = combined_data,
         mapping = aes(x = bio1, y = value, col = name))+
    geom_point()+
    geom_line()+
    facet_wrap(facets = "name",
               scales = "free",)
  
  ggplot(data = combined_data,
         mapping = aes(x = bio12, y = value, col = name))+
    geom_point()+
    geom_line()+
    facet_wrap(facets = "name",
               scales = "free",)
  
  
  ggplot(data = combined_data,
         mapping = aes(x = bio1, y = bio12, col = value))+
    geom_point()+
    facet_wrap(facets = "name",
               scales = "free",)
  
  p.list = lapply(sort(unique(combined_data$name)), function(i) {
    ggplot(data = combined_data[which(combined_data$name==i),],
           mapping = aes(x = bio1, y = bio12, col = value))+
      geom_point()+ggtitle(i)
    
    
      })
  
  library(gridExtra)
  do.call(grid.arrange, c(p.list, nrow=3))
  
    
  pp1 <- pp_gauss(p = pres_env$env[,1],
           bgrd = bg_env$env)

  #example plot  
  ex<-pp_gauss(p = as.data.frame(pres_env[,2]),
           bgrd = as.data.frame(bg_env[,2]))

  
  plot(y = dnorm(x = seq(from=10,to=40,by=0.1),
                 mean = ex$mean.p,
                 sd = ex$sigma.p)/
         dnorm(x = seq(from=10,to=40,by=0.1),
               mean = ex$mean.bgrd,
               sd = ex$sigma.bgrd),
       x = seq(from=10,to=40,by=0.1),
       pch=19,cex=0.1,
       col="purple",
       xlab="bio1",ylab="suitability (or density)",
       main="Gaussian")
  
  points(y = dnorm(x = seq(from=10,to=40,by=0.1),
                 mean = ex$mean.p,
                 sd = ex$sigma.p),
       x = seq(from=10,to=40,by=0.1),
       pch=19,cex=0.1,
       col="red")

  points(y = dnorm(x = seq(from=10,to=40,by=0.1),
                   mean = ex$mean.bgrd,
                   sd = ex$sigma.bgrd),
         x = seq(from=10,to=40,by=0.1),
         pch=19,cex=0.1,
         col="blue")
  abline(v = max(bg_env[,2],na.rm = T))
  abline(v = min(bg_env[,2],na.rm = T))
  abline(v = min(pres_env[,2],na.rm = T),lty=2)
  abline(v = max(pres_env[,2],na.rm = T),lty=2)
  
  #example kde
  
  ex2<-pp_kde(p = as.data.frame(pres_env[,2]),
               bgrd = as.data.frame(bg_env[,2]))
  
  
  plot(y = stats::fitted(np::npudens(bws=ex2$f1$bws,
                                     edat=data.frame(seq(from=10,to=40,by=0.01)),
                                     bwmethod='normal-reference'))/
         stats::fitted(np::npudens(bws=ex2$f0$bws,
                                   edat=data.frame(seq(from=10,to=40,by=0.01)),
                                   bwmethod='normal-reference')),
       x = seq(from=10,to=40,by=0.01),
       ylim=c(0,2),
       pch=19,cex=0.1,
       col="purple",
       xlab="bio1",ylab="suitability (or density)",main="KDE")
  
  points(y = stats::fitted(np::npudens(bws=ex2$f1$bws,
                                       edat=data.frame(seq(from=10,to=40,by=0.01)),
                                     bwmethod='normal-reference')),
         x = seq(from=10,to=40,by=0.01),
         pch=19,cex=.1,
         col="red")
  
  points(y = stats::fitted(np::npudens(bws=ex2$f0$bws,
                                     edat=data.frame(seq(from=10,to=40,by=0.01)),
                                     bwmethod='normal-reference')),
       x = seq(from=10,to=40,by=0.01),
       pch=19,cex=.1,
       col="blue")
  abline(v = max(bg_env[,2],na.rm = T))
  abline(v = min(bg_env[,2],na.rm = T))
  abline(v = min(pres_env[,2],na.rm = T),lty=2)
  abline(v = max(pres_env[,2],na.rm = T),lty=2)
  
  #example hybrid
  
  plot(y = dnorm(x = seq(from=10,to=40,by=0.01),
                 mean = ex$mean.p,
                 sd = ex$sigma.p)/
         stats::fitted(np::npudens(bws=ex2$f0$bws,
                                   edat=data.frame(seq(from=10,to=40,by=0.01)),
                                   bwmethod='normal-reference')),
       x = seq(from=10,to=40,by=0.01),
       pch=19,
       cex=0.1,
       col="purple",ylim=c(-0,6),
       xlab="bio1",ylab="suitability (or density)",main="hybrid")
  
  points(y = dnorm(x = seq(from=10,to=40,by=0.01),
                   mean = ex$mean.p,
                   sd = ex$sigma.p),
         x = seq(from=10,to=40,by=0.01),
         pch=19,
         cex=0.1,
         col="red")
  
  
  points(y = stats::fitted(np::npudens(bws=ex2$f0$bws,
                                       edat=data.frame(seq(from=10,to=40,by=0.01)),
                                       bwmethod='normal-reference')),
         x = seq(from=10,to=40,by=0.01),
         pch=19,cex=.1,
         col="blue")
  abline(v = max(bg_env[,2],na.rm = T))
  abline(v = min(bg_env[,2],na.rm = T))
  abline(v = min(pres_env[,2],na.rm = T),lty=2)
  abline(v = max(pres_env[,2],na.rm = T),lty=2)
  
  suit1<-
    PlugNPlay:::predict.pp_mod(object = pp1,
                               x = all_env_vals)
  
  pred <- env[[2]]
  pred <- setValues(pred,NA)
  
  pred[which(!na_or_not)] <- suit1
  plot(pred,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))
  
  
  pp2 <- PlugNPlay::pp_kde(p = pres_env,
                           bgrd = bg_env)
  
  suit2 <- PlugNPlay:::predict.pp_mod(object = pp2,
                                      x = all_env_vals) #interestingly, the predict here seems to be pretty slow.
  
  pred2 <- env[[2]]
  pred2 <- setValues(pred2,NA)
  
  pred2[which(!na_or_not)] <- suit2
  pred2[which(getValues(pred2)>10)] <- 0
  
  plot(pred2,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))
  pred2[which(getValues(pred2)>0.01)] <- 1
  plot(pred2)
  
#######################################################################################

  
#Gaussian

gauss_drake <-    
pp_gauss(p = pres_env,
         bgrd = bg_env,
         gauss_method = "classical")
  
gauss_suit_drake <-  
PlugNPlay:::predict.pp_mod(object = gauss_drake,
                           x = all_env_vals)  
  


gauss_me <-  
fit_plug_and_play(presence = pres_env$env,
                  background = bg_env$env,
                  method = "gaussian")

gauss_suit_me <- project_plug_and_play(pnp_model = gauss_me,
                                   data = bg_env$env)

gauss_pred_mine <- setValues(env[[2]],NA)
gauss_pred_mine[bg_env$bg_cells] <- gauss_suit_me
plot(gauss_pred_mine,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))

gauss_pred_drake <- env[[2]]
gauss_pred_drake <- setValues(gauss_pred_drake,NA)
gauss_pred_drake[which(!na_or_not)] <- gauss_suit_drake
plot(gauss_pred_drake,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))

###gauss checks out

# KDE
  #KDE projections are pretty slow, possibly due to dimensionality issue?

kde_drake <-    
  pp_kde(p = pres_env,
         bgrd = bg_env)
  
kde_suit_drake <-  
  PlugNPlay:::predict.pp_mod(object = kde_drake,
                             x = all_env_vals)  


kde_me <-  
  fit_plug_and_play(presence = pres_env,
                    background = bg_env,
                    method = "kde")

kde_suit_me <- project_plug_and_play(pnp_model = kde_me,
                                       data = all_env_vals)


kde_pred_mine <- env[[2]]
kde_pred_mine <- setValues(kde_pred_mine,NA)
kde_pred_mine[which(!na_or_not)] <- kde_suit_me
kde_pred_mine[which(getValues(kde_pred_mine)>10)]<-0
plot(kde_pred_mine,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))

kde_pred_drake <- env[[2]]
kde_pred_drake <- setValues(kde_pred_drake,NA)
kde_pred_drake[which(!na_or_not)] <- kde_suit_drake
kde_pred_drake[which(getValues(kde_pred_drake)>2)]<-0
plot(kde_pred_drake,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))

# kde checks out

####Hybrid

hybrid_me <-  
  fit_plug_and_play(presence = pres_env$env,
                    background = bg_env$env,
                    presence_method = "gaussian",
                    background_method = "gaussian",
                    bootstrap = "doublebag")

hybrid_suit_me <- project_plug_and_play(pnp_model = hybrid_me,
                                     data = all_env_vals)

gaussian_doublebag_auc <-evaluate_range_map(occurrences = occurrences,
                   env = env,
                   method = "gaussian",
                   bootstrap = "doublebag")

gaussian_auc <-evaluate_range_map(occurrences = occs[c("longitude","latitude")],
                                            env = env,
                                            method = "gaussian",
                                            bootstrap = "none",
                                  constraint_regions = continents)

gaussian_numbag <-evaluate_range_map(occurrences = occs[c("longitude","latitude")],
                                  env = env,
                                  method = "gaussian",
                                  bootstrap = "numbag",
                                  constraint_regions = continents)


#currently using np::npudens for kde, but there may be a faster method?
  #kdevine claims to not suffer from curse of dimensionality

hybrid_pred_mine <- setValues(env[[2]],NA)
hybrid_pred_mine[which(!na_or_not)] <- hybrid_suit_me
hybrid_pred_mine[which(getValues(hybrid_pred_mine) > 10)] <- 0
hybrid_pred_mine[which(getValues(hybrid_pred_mine) > 0.000001)] <- 1
plot(hybrid_pred_mine,xlim=c(-2000000,4000000),
     ylim=c(-1000000,3000000),
     main="hybrid")

plot(gauss_pred_mine,
     xlim=c(-2000000,4000000),
     ylim=c(-1000000,3000000),
     main="gaussian")

plot(kde_pred_mine,
     xlim=c(-2000000,4000000),
     ylim=c(-1000000,3000000),
     main="kde")


vine_model <- fit_plug_and_play(presence = pres_env,background = bg_env,method = "vine")

vine_suit <- project_plug_and_play(pnp_model = vine_model,
                                        data = all_env_vals)
vine_pred <- env[[2]]
vine_pred <- setValues(vine_pred,NA)
vine_pred[which(!na_or_not)] <- vine_suit
vine_pred[which(getValues(vine_pred) > 10)] <- 0
vine_pred[which(getValues(vine_pred) > 0.000001)] <- 1

plot(vine_pred)
plot(vine_pred,xlim=c(-2000000,4000000),
     ylim=c(-1000000,3000000),
     main="vine")
#######################################

num_hyb <-
fit_plug_and_play(presence = pres_env$env,
                  background = bg_env$env,
                  presence_method = "rangebagging",
                  background_method = "none")

pred_num_hyb <- project_plug_and_play(pnp_model = num_hyb,
                                      data = bg_env$env)

raster_pred_num_hyb <- setValues(env[[2]], NA)
raster_pred_num_hyb[which(!na_or_not)] <- pred_num_hyb
plot(raster_pred_num_hyb, main="rangebagging no numbag")
raster_pred_num_hyb[which(getValues(raster_pred_num_hyb) > 100)] <- 0
raster_pred_num_hyb[which(getValues(raster_pred_num_hyb) > 0.01)] <- 1

plot(raster_pred_num_hyb)
plot(raster_pred_num_hyb,
     xlim=c(-2000000,4000000),
     ylim=c(-1000000,3000000),
     main="vine")




#there are some copula methods that handle low-dimensional data
  #- cort is one, but doesn't seem to handle density estimation well.

#################################################
