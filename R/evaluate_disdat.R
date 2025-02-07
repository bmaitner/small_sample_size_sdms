#' @param presence_method See fit_plug_and_play
#' @param background_method See fit_plug_and_play
#' @param ratio_method see fit_density_ratio.  Only needed if presence and backgruond methods aren't supplied.
#' @param quantile Quantile for thresholding, set at 0.05
#' @param ncl Number of clusters to use for parallelizing. Defaults to 5, because 5-fold CV is what I parallelized
#' @importFrom pROC roc auc
#' @return List containing information on how well the selected model performs on the disdat datasets
evaluate_disdat <- function(presence_method = NULL,
                            background_method = NULL,
                            ratio_method = NULL,
                            quantile = 0.05,
                            verbose = TRUE,
                            ncl=5){
  
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
    
  if(verbose){message(paste("using presence/numerator method ",
                            presence_method, " and background/denominator method ",
                            background_method))}
    
    
  regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
  
  full_model_stats <- NULL
  fold_model_stats <- NULL
  
  cl <- makeCluster(ncl)
  registerDoParallel(cl)
  
  for (i in 1:length(regions)){
    
    
    region_i <- regions[i]
    data_i <- disData(region = region_i)
    
    if(verbose){message(paste("Starting region ",i, " of ", length(regions),": ",region_i))}
    

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
          
      if(verbose){message(paste("Starting species ",s, " of ",
                                length(unique(data_i$po$spid)),": ",species_s))
        }
      
      
      group_s <- unique(data_i$po$group[which(data_i$po$spid == species_s)])    
      presence_s <- data_i$po[which(data_i$po$spid == species_s),]
      background_s <- data_i$bg
      
      # re-scale presence and background.
      
        bg_means <- colMeans(background_s[,7:ncol(background_s)])
        bg_sd <- apply(X = background_s[,7:ncol(background_s)],MARGIN = 2,FUN = sd)
        
      

        presence_s[,7:ncol(presence_s)] <-
        S4DM:::rescale_w_objects(data = presence_s[,7:ncol(presence_s)],
                          mean_vector = bg_means,
                          sd_vector = bg_sd)
        
        background_s[,7:ncol(presence_s)] <-
          S4DM:::rescale_w_objects(data = background_s[,7:ncol(presence_s)],
                            mean_vector = bg_means,
                            sd_vector = bg_sd)


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
                        n_testing_background = NA,
                        runtime = NA,
                        entropy = NA)
      
      out_full <- data.frame(full_AUC = NA,
                             full_pAUC_specificity = NA,
                             full_pAUC_sensitivity = NA,
                             full_DOR = NA,
                             full_prediction_accuracy = NA,
                             full_sensitivity = NA,
                             full_specificity = NA,
                             full_correlation = NA,
                             full_kappa = NA,
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
                             runtime = NA,
                             entropy = NA)
      
      
      out <- foreach(fold = 1:length(unique(presence_data$fold)),
                     .packages = c("S4DM","tidyverse","DescTools"),
                     .combine = "rbind") %dopar% {
                       
                       if(verbose){message(paste("Starting fold ",fold, " of ",length(unique(presence_data$fold))))}
                       
                       # Populate fields that don't depend on model
                       
                       out$n_background[fold] <- nrow(background_s[,7:ncol(background_s)])
                       out$n_presence[fold] <- nrow(presence_s[which(presence_data$fold!=fold),
                                                               7:ncol(presence_s)])
                       out$n_testing_background[fold] <- nrow(background_s[,7:ncol(background_s)])
                       out$n_testing_presence[fold]  <- nrow(presence_s[which(presence_data$fold==fold),
                                                                        7:ncol(presence_s)])
                      
                       #This if statement skips cross validation if there is only one fold
                       
                       if(length(unique(presence_data$fold)) == 1){
                         print("Skipping cross validation, only one fold")
                         return(out)
                         
                       }
                       
                       model_fold <- NULL
                       
                       
                       if(is.null(ratio_method)){
                         
                         time_start <- Sys.time()              
                         
                         try(model_fold <- 
                               fit_plug_and_play(presence = presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)],
                                                 background = background_s[,7:ncol(background_s)],
                                                 presence_method = presence_method,
                                                 background_method = background_method),silent = T)
                         
                         time_finish <- Sys.time()
                         
                         model_time <- time_finish - time_start
                         
                         #convert model time to seconds if needed
                         
                         if(units(model_time) != "secs"){ units(model_time) <- "secs" }
                         
                         if(units(model_time) != "secs"){stop("Model time units not seconds")}
                         
                         model_time <- as.numeric(model_time)
                         
                       }else{
                         
                         time_start <- Sys.time()              
                         
                         
                         try(model_fold <- 
                               fit_density_ratio(presence = presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)],
                                                 background = background_s[,7:ncol(background_s)],
                                                 method = ratio_method),
                             silent = T)
                         
                         time_finish <- Sys.time()
                         
                         model_time <- time_finish - time_start
                         
                         #convert model time to seconds if needed
                         
                         if(units(model_time) != "secs"){ units(model_time) <- "secs" }
                         
                         if(units(model_time) != "secs"){stop("Model time units not seconds")}
                         
                         model_time <- as.numeric(model_time)
                         
                         
                       }
                       
                       
                       
                       
                       #if the fold model couldn't be fit, skip it (NA's will indicate this happened)
                       
                       if(is.null(model_fold)){
                         #next
                         
                         return(out %>%
                                  rename(fold_temp = fold) %>%
                                  filter(fold_temp == fold) %>%
                                  rename(fold = fold_temp)
                         ) #renaming to prevent issues with having the same variable name in environment and dataframe
                         
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
                       
                       
                       # Get entropy
                       
                       # setup needed files
                       fold_bg <- list()
                       fold_bg[[1]] <- background_s[,7:ncol(background_s)]
                       fold_bg[[2]] <- bg_means
                       fold_bg[[3]] <- bg_sd
                       names(fold_bg) <- c("env","env_mean","env_sd")
                       
                       fold_pres <- list()
                       fold_pres[[1]] <- presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)]
                       names(fold_pres) <- "env"
                       
                       
                       response_curves <- S4DM:::get_response_curves(env_bg = fold_bg,
                                                              env_pres = fold_pres,
                                                              pnp_model = model_fold,
                                                              n.int = 1000)
                       
                       mean_ent_fold <- 
                         response_curves %>%
                         select(variable,prediction) %>%
                         group_by(variable) %>%
                         summarise(entropy = Entropy(prediction)) %>%
                         ungroup() %>%
                         summarise(mean_ent = mean(entropy))
                       
                       class(mean_ent_fold)
                       
                       
                       fold_training_suitability_v_occurrence <- data.frame(suitability = training_predictions,
                                                                            occurrence = c(rep(1,length(which(presence_data$fold!=fold))),
                                                                                           rep(0,nrow(background_s))))
                       
                       fold_testing_suitability_v_occurrence <- data.frame(suitability = testing_predictions,
                                                                           occurrence = c(rep(1,length(which(presence_data$fold==fold))),
                                                                                          rep(0,nrow(background_s))))
                       
                       
                       training_roc_obj <- NULL
                       
                       training_roc_obj <- tryCatch(pROC::roc(response = fold_training_suitability_v_occurrence$occurrence,
                                                          predictor = fold_training_suitability_v_occurrence$suitability,
                                                          level = c(0,1),
                                                          direction = "<"),
                                                error = function(e){e})
                       
                       # only attempt to use ROC if it could be calculated properly
                       
                       if(inherits(training_roc_obj,"roc")){
                         
                         
                         out$training_AUC[fold] <- training_roc_obj$auc
                         
                         
                         out$training_pAUC_specificity[fold] <- pROC::auc(roc = training_roc_obj,
                                                                          partial.auc = c(.8, 1),
                                                                          partial.auc.correct = TRUE,
                                                                          partial.auc.focus = "specificity")[[1]]
                         
                         out$training_pAUC_sensitivity[fold] <- pROC::auc(roc = training_roc_obj,
                                                                          partial.auc = c(.8, 1),
                                                                          partial.auc.correct = TRUE,
                                                                          partial.auc.focus = "sensitivity")[[1]]
                         }
                         
                        
                       
                       #Testing data

                       testing_roc_obj <- NULL
                       
                       testing_roc_obj <- tryCatch(pROC::roc(response = fold_testing_suitability_v_occurrence$occurrence,
                                                              predictor = fold_testing_suitability_v_occurrence$suitability,
                                                              level = c(0,1),
                                                              direction = "<"),
                                                    error = function(e){e})
                       
                       
                       # only use ROC if it was correctly calculated
                       
                         if(inherits(testing_roc_obj,"roc")){
                           
                           out$testing_AUC[fold] <- testing_roc_obj$auc
                           
                           out$testing_pAUC_specificity[fold] <- pROC::auc(roc = testing_roc_obj,
                                                                           partial.auc = c(.8, 1),
                                                                           partial.auc.correct = TRUE,
                                                                           partial.auc.focus = "specificity")[[1]]
                           
                           out$testing_pAUC_sensitivity[fold] <- pROC::auc(roc = testing_roc_obj,
                                                                           partial.auc = c(.8, 1),
                                                                           partial.auc.correct = TRUE,
                                                                           partial.auc.focus = "sensitivity")[[1]]
                           
                           
                           
                         }
                         

                       
                       
                       # Code to make testing suitability scores binary
                       
                       threshold <- stats::quantile(x = fold_training_suitability_v_occurrence$suitability[which(fold_training_suitability_v_occurrence$occurrence==1)],
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
                       
                       out$runtime[fold] <- model_time
                       out$entropy[fold] <- mean_ent_fold$mean_ent
                       
                       return(out %>%
                                rename(fold_temp = fold) %>%
                                filter(fold_temp == fold) %>%
                                rename(fold = fold_temp)
                       ) #renaming to prevent issues with having the same variable name in environment and dataframe
                       
                     }#end fold
      
      
      
      #Fit full model  
      
      if(is.null(ratio_method)){
        
        time_start <- Sys.time()              
        
        model_full <- tryCatch(expr =  
                                 fit_plug_and_play(presence = presence_s[,7:ncol(presence_s)],
                                                   background = background_s[,7:ncol(background_s)],
                                                   presence_method = presence_method,
                                                   background_method = background_method),
                               error = function(e){
                                 return(NULL)
                               }
        )
        
        
        time_finish <- Sys.time()
        
        model_time_full <- time_finish - time_start
        
        if(units(model_time_full) != "secs"){ units(model_time_full) <- "secs" }
        
        if(units(model_time_full) != "secs"){stop("Model time full units not seconds")}
        
        model_time_full <- as.numeric(model_time_full)
        
        
        
        
      }else{
        
        time_start <- Sys.time()              
        
        model_full <- tryCatch(expr =  
                                 fit_density_ratio(presence = presence_s[,7:ncol(presence_s)],
                                                   background = background_s[,7:ncol(background_s)],
                                                   method = ratio_method),
                               error = function(e){
                                 return(NULL)
                               }
        )
        
        
        time_finish <- Sys.time()
        
        model_time_full <- time_finish - time_start
        
        if(units(model_time_full) != "secs"){ units(model_time_full) <- "secs" }
        
        if(units(model_time_full) != "secs"){stop("Model time full units not seconds")}
        
        model_time_full <- as.numeric(model_time_full)
        
        
        
        
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
      
      #get entropy
      
        # setup needed files
      
          full_bg <- list()
          full_bg[[1]] <- background_s[,7:ncol(background_s)]
          full_bg[[2]] <- bg_means
          full_bg[[3]] <- bg_sd
          names(full_bg) <- c("env","env_mean","env_sd")
        
          full_pres <- list()
          full_pres[[1]] <- presence_s[,7:ncol(presence_s)]
          names(full_pres) <- "env"
        
        
        response_curves_full <- S4DM:::get_response_curves(env_bg = full_bg,
                                               env_pres = full_pres,
                                               pnp_model = model_full,
                                               n.int = 1000)
      
          mean_ent_full <- 
            response_curves_full %>%
            select(variable,prediction) %>%
            group_by(variable) %>%
            summarise(entropy = Entropy(prediction)) %>%
            ungroup() %>%
            summarise(mean_ent = mean(entropy))
      
      
      full_suitability_v_occurrence <- 
        data.frame(suitability = full_predictions,
                   occurrence = c(rep(1,nrow(presence_s)),
                                  rep(0,nrow(background_s))))
      
      full_roc_obj <- NULL
      
      full_roc_obj <- tryCatch(pROC::roc(response = full_suitability_v_occurrence$occurrence,
                                               predictor = full_suitability_v_occurrence$suitability,
                                               level = c(0,1),
                                               direction = "<"),
               error = function(e){e})
      
      # only do auc stuff if the roc object worked properly
        
        if(inherits(full_roc_obj,"roc")){
          
        
        
        out_full$full_pAUC_specificity <- pROC::auc(roc = full_roc_obj,
                                              partial.auc = c(.8, 1),
                                              partial.auc.correct = TRUE,
                                              partial.auc.focus = "specificity")[[1]]
        
        out_full$full_pAUC_sensitivity <- pROC::auc(roc = full_roc_obj,
                                              partial.auc = c(.8, 1),
                                              partial.auc.correct = TRUE,
                                              partial.auc.focus = "sensitivity")[[1]]
        
        out_full$full_AUC <- full_roc_obj$auc
        
        }
      
      full_suitability_v_occurrence <- na.omit(full_suitability_v_occurrence)
      
      out_full$full_correlation <- cor(x = full_suitability_v_occurrence$suitability,
                                       y = full_suitability_v_occurrence$occurrence,
                                       method = "pearson")
      
      #Evaluate the full model with full model data
      
        threshold <- stats::quantile(x = full_suitability_v_occurrence$suitability[which(full_suitability_v_occurrence$occurrence == 1)],
                                     probs = quantile,
                                     na.rm = T)
      
        #Anything greater than suitability threshold is considered a presence  
          
          TP <- length(which(full_suitability_v_occurrence$suitability >= threshold &
                               full_suitability_v_occurrence$occurrence == 1))
          
          FN <- length(which(full_suitability_v_occurrence$suitability < threshold &
                               full_suitability_v_occurrence$occurrence == 1))
          
          TN <- length(which(full_suitability_v_occurrence$suitability < threshold &
                               full_suitability_v_occurrence$occurrence == 0))
          
          FP <- length(which(full_suitability_v_occurrence$suitability >= threshold &
                               full_suitability_v_occurrence$occurrence == 0))
          
        
        out_full$full_sensitivity <- TP / (TP + FN)
        out_full$full_specificity <- TN / (FP + TN)
        #precision <- TP / (TP + FP)
        out_full$full_DOR <- (TP*TN)/(FP*FN)
        #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
        out_full$full_prediction_accuracy <- (TP+TN)/(TP+TN+FP+FN)
        
        P_o <- (TP+TN)/(TP+TN+FP+FN)
        Ppres <- ((TP+FP)/(TP+TN+FP+FN))*((TP+FN)/(TP+TN+FP+FN))
        Pabs <- ((FN+TN)/(TP+TN+FP+FN))*((FP+TN)/(TP+TN+FP+FN))
        P_e <- Ppres+Pabs
        out_full$full_kappa <- (P_o - P_e)/(1-P_e)
      

      #Evaluate the full model with independent presence/absence data
      #For at least one location the PA data seem to be corrupted due to an error converting from wide to long format (the y column shows up as a species)
      
      
        pres_abs_data_s <- merge(x = data_i$pa[which(data_i$pa$spid == species_s),"siteid",drop=FALSE],
                                 y = data_i$env,
                                 sort = FALSE)
      
      # rescale presence abscence data
      
        pres_abs_data_s[,5:ncol(pres_abs_data_s)] <-
          S4DM:::rescale_w_objects(data = pres_abs_data_s[,5:ncol(pres_abs_data_s)],
                            mean_vector = bg_means,
                            sd_vector = bg_sd)
      
  
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
        
        pa_roc_obj <- NULL
        
        pa_roc_obj <- tryCatch( pROC::roc(response = pa_suitability_v_occurrence$occurrence,
                          predictor = pa_suitability_v_occurrence$suitability,
                          level = c(0,1),
                          direction = "<"),
                          error = function(e){e})
        
        # only do auc stuff if the roc object worked properly
        
        if(inherits(pa_roc_obj,"roc")){
          
        
        out_full$pa_pAUC_specificity <- pROC::auc(roc = pa_roc_obj,
                                            partial.auc = c(.8, 1),
                                            partial.auc.correct = TRUE,
                                            partial.auc.focus = "specificity")[[1]]
        
        out_full$pa_pAUC_sensitivity <- pROC::auc(roc = pa_roc_obj,
                                            partial.auc = c(.8, 1),
                                            partial.auc.correct = TRUE,
                                            partial.auc.focus = "sensitivity")[[1]]
        
        out_full$pa_AUC <- pa_roc_obj$auc
        
        }
        
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
        out_full$runtime <- model_time_full
        out_full$entropy <- mean_ent_full$mean_ent
      
      }#End code that is only run if the model was fit      
      
      #Save output
      
        full_model_stats <- rbind(full_model_stats,
                                  data.frame(species = species_s, out_full))
        fold_model_stats <- rbind(fold_model_stats,
                                  data.frame(species = species_s, out))
      

    }#s loop
    
    
  }#i loop
  
  # cleanup parallel stuff
  
  stopCluster(cl)
  rm(cl)
  
  # return output
  
  output <- list(full_model_stats = full_model_stats,
                 fold_model_stats = fold_model_stats)  
  
  return(output)    
  
  
}#end fx

