#' @param presence_method See fit_plug_and_play
#' @param background_method See fit_plug_and_play
#' @param ratio_method see fit_density_ratio.  Only needed if presence and backgruond methods aren't supplied.
#' @param quantile Quantile for thresholding, set at 0.05
#' @importFrom pROC roc auc
#' @return List containing information on how well the selected model performs on the disdat datasets
evaluate_disdat <- function(presence_method = NULL,
                            background_method = NULL,
                            ratio_method = NULL,
                            quantile = 0.05){
  
  if(is.null(presence_method) & is.null(background_method) & is.null(ratio_method)){
    stop("Please supply presence + background methods OR a ratio method")
   }
  
  if(!is.null(presence_method) & !is.null(background_method) & !is.null(ratio_method)){
    stop("Please ONLY supply presence + background methods OR a ratio method")
  }

  if(
    (!is.null(presence_method) & is.null(background_method)) | 
    (is.null(presence_method) & !is.null(background_method))
  ){
    
  stop("Please supply BOTH presence and background methods.")  
    
  }
    
  
  regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
  
  full_model_stats <- NULL
  fold_model_stats <- NULL
  
  for (i in 1:length(regions)){
    
    region_i <- regions[i]
    data_i <- disData(region = region_i)


    #Remove categorical predictors  
    
    if(region_i == "CAN"){
      data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("ontveg"))] 
      data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("ontveg"))]
      data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("ontveg"))]
      
    }
      
    if(region_i == "NSW"){
      data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("disturb","soilfert","vegsys"))] 
      data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("disturb","soilfert","vegsys"))]
      data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("disturb","soilfert","vegsys"))]
    }
    
    if(region_i == "NZ"){
      data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("age","toxicats"))] 
      data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("age","toxicats"))]
      data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("age","toxicats"))]
    }
    
    if(region_i == "SWI"){
      data_i$env <- data_i$env[which(!colnames(data_i$env) %in% c("calc","sfroyy"))] 
      data_i$bg <- data_i$bg[which(!colnames(data_i$bg) %in% c("calc","sfroyy"))]
      data_i$po <- data_i$po[which(!colnames(data_i$po) %in% c("calc","sfroyy"))]
    }
    
    
    #Get EPSG (way to not standardize...)
    if(region_i == "AWT"){epsg <- 28355}
    if(region_i == "CAN"){epsg <- 4008}
    if(region_i == "NSW"){epsg <- 4326}
    if(region_i == "NZ"){epsg <- 27200}
    if(region_i == "SA"){epsg <- 4326}
    if(region_i == "SWI"){epsg <- 21781}
    
    for(s in 1:length(unique(data_i$po$spid))){
      
      out <- NULL  
      
      species_s <- unique(data_i$po$spid)[s]
          
      group_s <- unique(data_i$po$group[which(data_i$po$spid == species_s)])    
      presence_s <- data_i$po[which(data_i$po$spid == species_s),]
      background_s <- data_i$bg
      
      #stratify presences
      
        presence_data <-
          stratify_spatial(occurrence_sf = st_as_sf(x = presence_s[c("x","y")],
                                                    coords = c("x","y")) |>
                             st_set_crs(epsg),
                           nfolds = NULL,
                           nsubclusters = NULL)

      #Make empty output
      
      out <- data.frame(fold = 1:length(unique(presence_data$fold)),
                        training_AUC = NA,
                        training_pAUC_specificity = NA,
                        training_pAUC_sensitivity = NA,
                        testing_AUC = NA,
                        testing_pAUC_specificity = NA,
                        testing_pAUC_sensitivity = NA,
                        testing_DOR = NA,
                        testing_prediction_accuracy = NA,
                        testing_sensitivity = NA,
                        testing_specificity = NA,
                        testing_correlation = NA,
                        testing_kappa = NA,
                        n_presence = NA,
                        n_background = NA,
                        n_testing_presence = NA,
                        n_testing_background = NA)
      
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
                             n_pa_absence = NA)
      
      
      for(fold in 1:length(unique(presence_data$fold))){
        
        #This if statement skips cross validation if there is only one fold
        if(length(unique(presence_data$fold)) == 1){
          print("Skipping cross validation, only one fold")
          next
          }
        
        
        model_fold <- NULL
        
        
        if(is.null(ratio_method)){
          
          try(model_fold <- 
                fit_plug_and_play(presence = presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)],
                                  background = background_s[,7:ncol(background_s)],
                                  presence_method = presence_method,
                                  background_method = background_method),silent = T)  
          
        }else{
          
          try(model_fold <- 
                fit_density_ratio(presence = presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)],
                                  background = background_s[,7:ncol(background_s)],
                                  method = ratio_method),
              silent = T)
          
        }
        
        
        
        
        #if the fold model couldn't be fit, skip it (NA's will indicate this happened)
        
        if(is.null(model_fold)){
          next
        }
        
        
        pres_vector <- paste(presence_s$x[which(presence_data$fold!=fold)],
                             presence_s$y[which(presence_data$fold!=fold)])
        
        bg_vector <- paste(background_s$x,
                           background_s$y) 
        
        if(any(pres_vector %in% bg_vector)){stop("write more code")}
        
        training_data <- rbind(presence_s[which(presence_data$fold!=fold),
                                          7:ncol(presence_s)],
                               background_s[,7:ncol(background_s)])
        
        testing_data <- rbind(presence_s[which(presence_data$fold==fold),
                                         7:ncol(presence_s)],
                              background_s[,7:ncol(background_s)])
        
        
        if(is.null(ratio_method)){
          
          training_predictions <- project_plug_and_play(pnp_model = model_fold,
                                                        data = training_data)
          
          testing_predictions <- project_plug_and_play(pnp_model = model_fold,
                                                       data = testing_data)
          
        }else{
            
          training_predictions <- project_density_ratio(dr_model = model_fold,
                                                        data = training_data)
          
          testing_predictions <- project_density_ratio(dr_model = model_fold,
                                                       data = testing_data)
          
          
          }
        
        
        
        fold_training_suitability_v_occurrence <- data.frame(suitability = training_predictions,
                                                             occurrence = c(rep(1,length(which(presence_data$fold!=fold))),
                                                                            rep(0,nrow(background_s))))
        
        fold_testing_suitability_v_occurrence <- data.frame(suitability = testing_predictions,
                                                            occurrence = c(rep(1,length(which(presence_data$fold==fold))),
                                                                           rep(0,nrow(background_s))))
        
        
        training_roc_obj <- pROC::roc(response = fold_training_suitability_v_occurrence$occurrence,
                                predictor = fold_training_suitability_v_occurrence$suitability)
        
        out$training_AUC[fold] <- training_roc_obj$auc
        
        
        
        out$training_pAUC_specificity[fold] <- pROC::auc(roc = training_roc_obj,
                                                partial.auc = c(.8, 1),
                                                partial.auc.correct = TRUE,
                                                partial.auc.focus = "specificity")[[1]]
        
        out$training_pAUC_sensitivity[fold] <- pROC::auc(roc = training_roc_obj,
                                                partial.auc = c(.8, 1),
                                                partial.auc.correct = TRUE,
                                                partial.auc.focus = "sensitivity")[[1]]
        
        #Testing data
        
        testing_roc_obj <- pROC::roc(response = fold_testing_suitability_v_occurrence$occurrence,
                               predictor = fold_testing_suitability_v_occurrence$suitability)
        
        out$testing_AUC[fold] <- testing_roc_obj$auc
        
        
        
        out$testing_pAUC_specificity[fold] <- pROC::auc(roc = testing_roc_obj,
                                               partial.auc = c(.8, 1),
                                               partial.auc.correct = TRUE,
                                               partial.auc.focus = "specificity")[[1]]
        
        out$testing_pAUC_sensitivity[fold] <- pROC::auc(roc = testing_roc_obj,
                                               partial.auc = c(.8, 1),
                                               partial.auc.correct = TRUE,
                                               partial.auc.focus = "sensitivity")[[1]]
        
        
        # Code to make testing suitability scores binary
        
        threshold <- stats::quantile(x = fold_testing_suitability_v_occurrence$suitability[which(fold_testing_suitability_v_occurrence$occurrence==1)],
                                     probs = quantile,
                                     na.rm = T)
        
        #Anything greater than suitability threshold is considered a presence  
        
        TP <- length(which(fold_testing_suitability_v_occurrence$suitability >= threshold &
                             fold_testing_suitability_v_occurrence$occurrence == 1))
        
        FN <- length(which(fold_testing_suitability_v_occurrence$suitability < threshold &
                             fold_testing_suitability_v_occurrence$occurrence == 1))
        
        TN <- length(which(fold_testing_suitability_v_occurrence$suitability < threshold &
                             fold_testing_suitability_v_occurrence$occurrence == 0))
        
        FP <- length(which(fold_testing_suitability_v_occurrence$suitability >= threshold &
                             fold_testing_suitability_v_occurrence$occurrence == 0))
        
        
        
        
        sensitivity <- TP / (TP + FN)
        specificity <- TN / (FP + TN)
        #precision <- TP / (TP + FP)
        DOR <- (TP*TN)/(FP*FN)
        #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
        prediction_accuracy <- (TP+TN)/(TP+TN+FP+FN)
        P_o <- (TP+TN)/(TP+TN+FP+FN)
        Ppres <- ((TP+FP)/(TP+TN+FP+FN))*((TP+FN)/(TP+TN+FP+FN))
        Pabs <- ((FN+TN)/(TP+TN+FP+FN))*((FP+TN)/(TP+TN+FP+FN))
        P_e <- Ppres+Pabs
        kappa <- (P_o - P_e)/(1-P_e)
        
        
        out$testing_DOR[fold] <- DOR
        out$testing_prediction_accuracy[fold] <- prediction_accuracy
        out$testing_sensitivity[fold] <- sensitivity
        out$testing_specificity[fold] <- specificity
        out$testing_kappa[fold] <- kappa
        
        fold_testing_suitability_v_occurrence <- na.omit(fold_testing_suitability_v_occurrence)
        out$testing_correlation[fold] <- cor(fold_testing_suitability_v_occurrence$suitability,fold_testing_suitability_v_occurrence$occurrence)
        
        out$n_background[fold] <- nrow(background_s[,7:ncol(background_s)])
        out$n_presence[fold] <- nrow(presence_s[which(presence_data$fold!=fold),
                                  7:ncol(presence_s)])
        out$n_testing_background[fold] <- nrow(background_s[,7:ncol(background_s)])
        out$n_testing_presence[fold]  <- nrow(presence_s[which(presence_data$fold==fold),
                                                   7:ncol(presence_s)])
        
        
      }#end fold
      
      
      
      #Fit full model  
      
      if(is.null(ratio_method)){
        
        model_full <- tryCatch(expr =  
                                 fit_plug_and_play(presence = presence_s[,7:ncol(presence_s)],
                                                   background = background_s[,7:ncol(background_s)],
                                                   presence_method = presence_method,
                                                   background_method = background_method),
                               error = function(e){
                                 return(NULL)
                               }
        )        
        
        
        
      }else{
        
        
        model_full <- tryCatch(expr =  
                                 fit_density_ratio(presence = presence_s[,7:ncol(presence_s)],
                                                   background = background_s[,7:ncol(background_s)],
                                                   method = ratio_method),
                               error = function(e){
                                 return(NULL)
                               }
        )
        
        
        
        
        
      }
      

      
      
      
      
      
      
      
      
      
      
      
      if(!is.null(model_full)){
      
      full_data <- rbind(presence_s[,7:ncol(presence_s)],
                         background_s[,7:ncol(background_s)])
      
      
      if(is.null(ratio_method)){
        
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
      
      
      pres_abs_data_s <- merge(x = data_i$pa[which(data_i$pa$spid == species_s),"siteid",drop=FALSE],
                               y = data_i$env,
                               sort = FALSE)
      
      if(!all(pres_abs_data_s$siteid == data_i$pa$siteid[which(data_i$pa$spid == species_s)])){
        stop("Problem with data order in P/A data")
      }
      

      #Estimate suitabilities for PA data
      
      
      if(is.null(ratio_method)){
        
        pa_predictions <- project_plug_and_play(pnp_model = model_full,
                                                data = pres_abs_data_s[,5:ncol(pres_abs_data_s)])  
        
      }else{
        
        pa_predictions <- project_density_ratio(dr_model = model_full,
                                                data = pres_abs_data_s[,5:ncol(pres_abs_data_s)])
        
      }
      
      
      pa_suitability_v_occurrence <- 
        data.frame(suitability = pa_predictions,
                   occurrence =  data_i$pa$pa[which(data_i$pa$spid == species_s)])
      
      
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
      
      }#End code that is only run if the model was fit      
      
      #Save output
      full_model_stats <- rbind(full_model_stats,
                                data.frame(species = species_s, out_full))
      fold_model_stats <- rbind(fold_model_stats,
                                data.frame(species = species_s, out))
      
      
      
    }#s loop
    
    
  }#i loop
  
  
  output <- list(full_model_stats = full_model_stats,
                 fold_model_stats = fold_model_stats)  
  
  return(output)    
  
  
}#end fx

