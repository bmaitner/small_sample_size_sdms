

library(disdat)
library(tsutils)
data_i <- disData(region = "NSW")
data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("disturb","soilfert","vegsys"))] 
data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("disturb","soilfert","vegsys"))]
data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("disturb","soilfert","vegsys"))]
epsg <- 4326


species_s <- unique(data_i$po$spid)[9]
group_s <- unique(data_i$po$group[which(data_i$po$spid == species_s)])    
presence_s <- data_i$po[which(data_i$po$spid == species_s),]
background_s <- data_i$bg

p.f <- c(rep(1,nrow(presence_s)),rep(0,nrow(background_s)))
data.f <- rbind(presence_s, background_s)[7:ncol(presence_s)]

data.f <- data.f[1:1000,]
p.f <- p.f[1:1000]  


 model <- pbsdm::fit_plug_and_play(presence = presence_s[7:ncol(presence_s)],
                                   background = background_s[7:ncol(presence_s)],
                                   method = "gaussian")

 
 
 predictions <- project_plug_and_play(pnp_model = model,data = data.f)

 binary_preds <- threshold_predictions(predictions = predictions,
                       presences = p.f)


 
 ##########################################
 
 presences <- presence_s[7:ncol(presence_s)]
 background <- background_s[7:ncol(presence_s)]
 
 fit_plug_and_play()

 
#' @param presences covariates at presence locations
#' @param background covariates at backgruond locations
#' @param selection_criterion to use for ranking
#' @param return_model_info If TRUE, returns a data.frame of model fit info instead of just the best option
tune_gaussian <- function(presences,
                           background,
                           selection_criterion = "prediction_accuracy",
                           return_model_info = F){
  
   #' @param type one of either "classical", "robust", or "regularized" (the default)
   classical_fit <- evaluate_model(presences = presences,
                                   background = background,
                                   method = "gaussian",
                                   type = "classical")

   robust_fit <- evaluate_model(presences = presences,
                                   background = background,
                                   method = "gaussian",
                                   type = "robust")
   
   regularized_fit <- evaluate_model(presences = presences,
                                background = background,
                                method = "gaussian",
                                type = "regularized")
   
   out <- rbind(classical_fit,robust_fit,regularized_fit)
   out$type <- c("classical","robust","regularized")
   colnames(out) <- gsub(pattern = "full_",replacement = "",x = colnames(out))
   
   if(!selection_criterion %in% colnames(out)){
     message("selection_criterion doesn't match returned column names, returning NA")
     
     return(data.frame(type = NA))
   }
   
   
   if(!return_model_info){
     
     return(data.frame(type = out$type[which.max(out[,selection_criterion])]))
    
   }
   
   return(out)
   
   
   
 }


#############################
