

# get a set of species occurrence data

  occurrence_files <-
  list.files(path = "C:/Users/Brian Maitner/Desktop/current_projects/BIENWorkflow/inst/extdata/Demo_SpeciesCSVs/",
             pattern = ".csv",
             full.names = T)

# get a set of climate data
  library(raster)
  env <- stack(x = "C:/Users/Brian Maitner/Desktop/current_projects/BIENWorkflow/inst/extdata/Demo_Env/AllEnv.tif")
  names(env) <- read.csv("C:/Users/Brian Maitner/Desktop/current_projects/BIENWorkflow/inst/extdata/Demo_Env/layerNames.csv",header = F)[,1]

# load Drake Lab's plug n play package to check out
  library(PlugNPlay)
  
  
##################

  #Quick and dirty sampling of presences and absences

  source("R/get_env_pres.R")
  source("R/get_env_bg.R")
  source("R/fit_plug_and_play.R")
  source("R/pnp_gaussian.R")
  source("R/project_plug_and_play.R")
  
  all_env <- getValues(env)
  all_env_vals <-na.omit(all_env)
  
  na_or_not<-
    apply(X = all_env,
          MARGIN = 1,
          FUN = function(x){
            any(is.na(x))
            
          }
    )
  
  occs <- read.csv(occurrence_files[3])
  
  pres_env <- na.omit(get_env_pres(coords = occs[c("longitude","latitude")],
                           env = env)) 
  bg_env <- na.omit(get_env_bg(coords = occs[c("longitude","latitude")],
                       env = env,
                       width = 50000)) 
  
  pp1 <- pp_gauss(p = pres_env[,1],
           bgrd = bg_env)

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
  
  ###################################
    suit1<-
    PlugNPlay:::predict.pp_mod(object = pp1,
                             x = all_env_vals)
  
  pred <- env[[2]]
  pred <- setValues(pred,NA)
  
  pred[which(!na_or_not)] <- suit1
  plot(pred,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))
    
hist(all_env_vals[,4])  
  
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

#make example plot to wrap head around idea
  
  
#######################################################################################
#######################################################################################

  



test_model <-  
fit_plug_and_play(presence = pres_env,
                  background = bg_env,
                  method = "gaussian")

test_suit <- project_plug_and_play(pnp_model = test_model,
                                   data = all_env_vals)

pred_mine <- env[[2]]
pred_mine <- setValues(pred_mine,NA)

pred_mine[which(!na_or_not)] <- test_suit
plot(pred,xlim=c(-2000000,4000000),ylim=c(-1000000,3000000))
