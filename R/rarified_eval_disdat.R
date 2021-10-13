# Need to figure out a good method for rarifying points

#start by histograms of number of presence points

presences <- unique(full_model_outputs[c("species","n_presence")])

(2:10)^2

stop("Coding still in progress")

model_vector <- unique(full_model_output_all$method)


library(disdat)
library(tidyverse)
library(pbsdm)
#' @param presence_vector Vector of presence values to rarify things to.  All species with values greater than or equal to the max value will be included.
#' @param n_reps Number of replicates for each presence level
#' @param model_vector A vector of model names.  For presence/background, these should be specified as presence/background.
#' @param quantile Quantile for thresholding, set at 0.05
rarified_eval_disdat <- function(presence_vector = (2:20)^2,
                                 n_reps = 3,
                                 model_vector,
                                 quantile = 0.05,
                                 temp_RDS = "outputs/temp_rarified.RDS"){

  full_model_stats <- NULL
  # First, get a list of relevant species
  
        regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
        
        presence_counts <- NULL
        
        for (i in 1:length(regions)){
          
          region <- regions[i]
          data_i <- disData(region = region)  
          
      
          data_i$po %>%
            group_by(spid) %>%
            summarize(po_presence = n()) %>%
            mutate(region = region) %>% 
            rbind(presence_counts) -> presence_counts
      
          
        }
  
  # Get the subset that meet the required criteria
        
      presence_counts %>%
        filter(po_presence >= max(presence_vector)) -> focal_species
      
      cat("Using", nrow(focal_species), "of ", nrow(presence_counts), "species that meet the presence number criteria")
      
      
  # Iterate through regions,species,reps, and models.    
      
  
    for(r in 1:length(unique(focal_species$region))){
      
          region <- unique(focal_species$region)[r]
      
          data_i <- disData(region = region)
          
          
          #Remove categorical predictors  
          
          if(region == "CAN"){
            data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("ontveg"))] 
            data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("ontveg"))]
            data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("ontveg"))]
            
          }
          
          if(region == "NSW"){
            data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("disturb","soilfert","vegsys"))] 
            data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("disturb","soilfert","vegsys"))]
            data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("disturb","soilfert","vegsys"))]
          }
          
          if(region == "NZ"){
            data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("age","toxicats"))] 
            data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("age","toxicats"))]
            data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("age","toxicats"))]
          }
          
          if(region == "SWI"){
            data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("calc","sfroyy"))] 
            data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("calc","sfroyy"))]
            data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("calc","sfroyy"))]
          }
          
          
          #Get EPSG (way to not standardize...)
          if(region == "AWT"){epsg <- 28355}
          if(region == "CAN"){epsg <- 4008}
          if(region == "NSW"){epsg <- 4326}
          if(region == "NZ"){epsg <- 27200}
          if(region == "SA"){epsg <- 4326}
          if(region == "SWI"){epsg <- 21781}
          
          
      for(s in 1:length(unique(focal_species$spid))){
        
        
            #Pull the necessary species data
            species <- unique(focal_species$spid)[s]
            group_s <- unique(data_i$po$group[which(data_i$po$spid == species)])    
            #presence_s <- data_i$po[which(data_i$po$spid == species),]
            background_s <- data_i$bg

        for(p in presence_vector){
          for(n in 1:n_reps){
            
            presence_s <- data_i$po[which(data_i$po$spid == species),]
            presence_s <- presence_s[sample(x = nrow(presence_s), size = p, replace = F),]
            
            
            for(m in 1:length(unique(model_vector))){
  
              
              
              
              replicate <- n
              model <- unique(model_vector)[m]
            
              
              #Parse the model text
                pres_method <- NULL
                bg_method <- NULL
                dr_method <- NULL
                
                #If a p-b model
                if(grepl(pattern = "/",x = model)){
                  
                  model <- gsub(pattern = " ",replacement = "",x = model)
                  pres_method <- strsplit(x = model,split = "/")[[1]][1]
                  bg_method <- strsplit(x = model,split = "/")[[1]][2]
                  
                }else{
                
                  dr_method <- model  
                  
                }
              
              
                #Fit the model
                
                #Fit full model  
                
                if(is.null(dr_method)){
                  model_full <- NULL
                  runtime <- NULL
                  runtime <- proc.time() 
                  
                  
                  
                  model_full <- tryCatch(expr =  
                                           fit_plug_and_play(presence = presence_s[,7:ncol(presence_s)],
                                                             background = background_s[,7:ncol(background_s)],
                                                             presence_method = pres_method,
                                                             background_method = bg_method),
                                         error = function(e){
                                           message("problem fitting")
                                           return(NULL)
                                         }
                  )        
                  runtime <- proc.time() - runtime
                  
                  
                }else{
                  
                  model_full <- NULL
                  runtime <- NULL
                  
                  runtime <- proc.time() 
                  
                  model_full <- tryCatch(expr =  
                                           fit_density_ratio(presence = presence_s[,7:ncol(presence_s)],
                                                             background = background_s[,7:ncol(background_s)],
                                                             method = dr_method),
                                         error = function(e){
                                           return(NULL)
                                         }
                  )
                  runtime <- proc.time() - runtime
                  
                }
                
                
                out_full <- data.frame(full_AUC = NA,
                                       full_pAUC_specificity = NA,
                                       full_pAUC_sensitivity = NA,
                                       full_correlation = NA,
                                       pa_AUC = NA,
                                       pa_pAUC_specificity = NA,
                                       pa_pAUC_sensitivity = NA,
                                       pa_DOR = NA,
                                       pa_prediction_accuracy = NA,
                                       pa_sensitivity = NA,
                                       pa_specificity = NA,
                                       pa_correlation = NA,
                                       pa_kappa = NA,
                                       n_presence = NA,
                                       n_background = NA,
                                       n_pa_presence = NA,
                                       n_pa_absence = NA,
                                       runtime = NA)
                
                
                #If the model failed, just record NAs and metadata
                
            
                if(is.null(model_full)){
                  
                  out_full$n_background <- nrow(background_s[,7:ncol(background_s)])
                  out_full$n_presence <- nrow(presence_s[,7:ncol(presence_s)])
                  out_full$n_pa_absence <- length(data_i$pa$pa[which(data_i$pa$spid == species & data_i$pa$pa == 0)])
                  out_full$n_pa_presence <- length(data_i$pa$pa[which(data_i$pa$spid == species & data_i$pa$pa == 1)])
                  out_full$runtime <- runtime[3]
                  
                  
                  full_model_stats <- rbind(full_model_stats,
                                            data.frame(species = species,
                                                       model = model, out_full))
                  
                  saveRDS(object = full_model_stats,file = temp_RDS)
                  
                 next 
                  
                }
                
                
                if(!is.null(model_full)){
                  
                  full_data <- rbind(presence_s[,7:ncol(presence_s)],
                                     background_s[,7:ncol(background_s)])
                  
                  
                  if(is.null(dr_method)){
                    
                    full_predictions <- project_plug_and_play(pnp_model = model_full,
                                                              data = full_data)  
                    
                  }else{
                    
                    full_predictions <- project_density_ratio(dr_model = model_full,
                                                              data = full_data)
                    
                  }
                  
                  
                  
                  full_suitability_v_occurrence <- 
                    data.frame(suitability = full_predictions,
                               occurrence = c(rep(1,nrow(presence_s)),
                                              rep(0,nrow(background_s))))
                  
                  full_roc_obj <- pROC::roc(response = full_suitability_v_occurrence$occurrence,
                                            predictor = full_suitability_v_occurrence$suitability)
                  
                  
                  
                  out_full$full_pAUC_specificity <- pROC::auc(roc = full_roc_obj,
                                                              partial.auc = c(.8, 1),
                                                              partial.auc.correct = TRUE,
                                                              partial.auc.focus = "specificity")[[1]]
                  
                  out_full$full_pAUC_sensitivity <- pROC::auc(roc = full_roc_obj,
                                                              partial.auc = c(.8, 1),
                                                              partial.auc.correct = TRUE,
                                                              partial.auc.focus = "sensitivity")[[1]]
                  
                  out_full$full_AUC <- full_roc_obj$auc
                  
                  full_suitability_v_occurrence <- na.omit(full_suitability_v_occurrence)
                  
                  out_full$full_correlation <- cor(x = full_suitability_v_occurrence$suitability,
                                                   y = full_suitability_v_occurrence$occurrence,
                                                   method = "pearson")
                  
                  
                  #Evaluate the full model with independent presence/absence data
                  #For at least one location the PA data seem to be corrupted due to an error converting from wide to long format (the y column shows up as a species)
                  
                  
                  pres_abs_data_s <- merge(x = data_i$pa[which(data_i$pa$spid == species),"siteid",drop=FALSE],
                                           y = data_i$env,
                                           sort = FALSE)
                  
                  if(!all(pres_abs_data_s$siteid == data_i$pa$siteid[which(data_i$pa$spid == species)])){
                    stop("Problem with data order in P/A data")
                  }
                  
                  
                  #Estimate suitabilities for PA data
                  
                  
                  if(is.null(dr_method)){
                    
                    pa_predictions <- project_plug_and_play(pnp_model = model_full,
                                                            data = pres_abs_data_s[,5:ncol(pres_abs_data_s)])  
                    
                  }else{
                    
                    pa_predictions <- project_density_ratio(dr_model = model_full,
                                                            data = pres_abs_data_s[,5:ncol(pres_abs_data_s)])
                    
                  }
                  
                  
                  pa_suitability_v_occurrence <- 
                    data.frame(suitability = pa_predictions,
                               occurrence =  data_i$pa$pa[which(data_i$pa$spid == species)])
                  
                  
                  pa_roc_obj <- pROC::roc(response = pa_suitability_v_occurrence$occurrence,
                                          predictor = pa_suitability_v_occurrence$suitability)
                  
                  out_full$pa_pAUC_specificity <- pROC::auc(roc = pa_roc_obj,
                                                            partial.auc = c(.8, 1),
                                                            partial.auc.correct = TRUE,
                                                            partial.auc.focus = "specificity")[[1]]
                  
                  out_full$pa_pAUC_sensitivity <- pROC::auc(roc = pa_roc_obj,
                                                            partial.auc = c(.8, 1),
                                                            partial.auc.correct = TRUE,
                                                            partial.auc.focus = "sensitivity")[[1]]
                  
                  out_full$pa_AUC <- pa_roc_obj$auc
                  
                  pa_suitability_v_occurrence <- na.omit(pa_suitability_v_occurrence)
                  
                  out_full$pa_correlation <- cor(x = pa_suitability_v_occurrence$suitability,
                                                 y = pa_suitability_v_occurrence$occurrence,
                                                 method = "pearson")
                  
                  # Code to make testing suitability scores binary
                  
                  #For this one, we set the threshold based on the full model using presence and background
                  
                  
                  threshold <- stats::quantile(x = full_suitability_v_occurrence$suitability[which(full_suitability_v_occurrence$occurrence==1)],
                                               probs = quantile,
                                               na.rm = T)
                  
                  #Anything greater than suitability threshold is considered a presence  
                  
                  TP <- length(which(pa_suitability_v_occurrence$suitability >= threshold &
                                       pa_suitability_v_occurrence$occurrence == 1))
                  
                  FN <- length(which(pa_suitability_v_occurrence$suitability < threshold &
                                       pa_suitability_v_occurrence$occurrence == 1))
                  
                  TN <- length(which(pa_suitability_v_occurrence$suitability < threshold &
                                       pa_suitability_v_occurrence$occurrence == 0))
                  
                  FP <- length(which(pa_suitability_v_occurrence$suitability >= threshold &
                                       pa_suitability_v_occurrence$occurrence == 0))
                  
                  
                  
                  
                  out_full$pa_sensitivity <- TP / (TP + FN)
                  out_full$pa_specificity <- TN / (FP + TN)
                  #precision <- TP / (TP + FP)
                  out_full$pa_DOR <- (TP*TN)/(FP*FN)
                  #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
                  out_full$pa_prediction_accuracy <- (TP+TN)/(TP+TN+FP+FN)
                  
                  P_o <- (TP+TN)/(TP+TN+FP+FN)
                  Ppres <- ((TP+FP)/(TP+TN+FP+FN))*((TP+FN)/(TP+TN+FP+FN))
                  Pabs <- ((FN+TN)/(TP+TN+FP+FN))*((FP+TN)/(TP+TN+FP+FN))
                  P_e <- Ppres+Pabs
                  out_full$pa_kappa <- (P_o - P_e)/(1-P_e)
                  
                  out_full$n_background <- nrow(background_s[,7:ncol(background_s)])
                  out_full$n_presence <- nrow(presence_s[,7:ncol(presence_s)])
                  out_full$n_pa_absence <- length(which(pa_suitability_v_occurrence$occurrence == 0))
                  out_full$n_pa_presence <- length(which(pa_suitability_v_occurrence$occurrence == 1))
                  out_full$runtime <- runtime[3]
                  
                }#End code that is only run if the model was fit      
                
                #Save output
                full_model_stats <- rbind(full_model_stats,
                                          data.frame(species = species,
                                                     model = model, out_full))
                
                saveRDS(object = full_model_stats,file = temp_RDS)
                
              
              
            } #m loop
          } #n loop
        } #p presences loop
      }#s loop
    } #r loop 
      
      
  
  
  
  
  
}# end function














