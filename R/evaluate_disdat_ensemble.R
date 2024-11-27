#' @param presence_method See fit_plug_and_play
#' @param background_method See fit_plug_and_play
#' @param ratio_method see fit_density_ratio.  Only needed if presence and backgruond methods aren't supplied.
#' @param quantile Quantile for thresholding, set at 0.05
#' @param ncl Number of clusters to use for parallelizing. Defaults to 5, because 5-fold CV is what I parallelized
#' @importFrom pROC roc auc
#' @return List containing information on how well the selected model performs on the disdat datasets
evaluate_ensemble_disdat <- function(model_vector = NULL,
                            quantile = 0.05,
                            verbose = TRUE,
                            ncl=5,
                            temp_file = "outputs/temp_ensemble_eval.RDS"){

  regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
  
  full_model_stats <- NULL
  fold_model_stats <- NULL
  
  cl <- makeCluster(ncl)
  registerDoParallel(cl)
  
  if(file.exists(temp_file)){
    
    full_model_stats <- readRDS(temp_file)$full_model_stats
    fold_model_stats <- readRDS(temp_file)$fold_model_stats
    
    
    }
  
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
      
    # Skip the species if its already been done
      
      if(species_s %in% full_model_stats$species &
         species_s %in% fold_model_stats$species){
        next
        }
          
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
        pbsdm:::rescale_w_objects(data = presence_s[,7:ncol(presence_s)],
                          mean_vector = bg_means,
                          sd_vector = bg_sd)
        
        background_s[,7:ncol(presence_s)] <-
          pbsdm:::rescale_w_objects(data = background_s[,7:ncol(presence_s)],
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
                          training_vote_AUC = NA,
                          
                          training_pAUC_specificity = NA,
                          training_vote_pAUC_specificity = NA,
                          
                          training_pAUC_sensitivity = NA,
                          training_vote_pAUC_sensitivity = NA,
                          
                          testing_AUC = NA,
                          testing_vote_AUC = NA,
                          
                          testing_pAUC_specificity = NA,
                          testing_vote_pAUC_specificity = NA,
                          
                          testing_pAUC_sensitivity = NA,
                          testing_vote_pAUC_sensitivity = NA,
                          
                          testing_DOR = NA,
                          testing_vote_DOR = NA,
                          
                          testing_prediction_accuracy = NA,
                          testing_vote_prediction_accuracy = NA,
                          
                          testing_sensitivity = NA,
                          testing_vote_sensitivity = NA,
                          
                          testing_specificity = NA,
                          testing_vote_specificity = NA,
                          
                          testing_correlation = NA,
                          testing_vote_correlation = NA,
                          
                          testing_kappa = NA,
                          testing_vote_kappa = NA,
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
      
      
      if(verbose){message("Starting CV")}
      
      out <- foreach(fold = 1:length(unique(presence_data$fold)),
                     .packages = c("pbsdm","tidyverse","DescTools"),
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
                         message("Skipping cross validation, only one fold")
                         return(out)
                         
                       }
                       

                    ### ensemble fitting
                       time_start <- Sys.time()              
                       
                       model_fold <- list()
                       
                       
                       for(m in 1:length(model_vector)){
                         
                         model_m <- model_vector[m] 
                         
                         model_m <- gsub(pattern = " ",replacement = "",x = model_m)
                         
                         split_m <- strsplit(x = model_m,split = "/") %>% unlist()
                         
                         if(length(split_m) == 1){
                           
                           
                           try(model_fold[[m]] <- 
                                 fit_density_ratio(presence = presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)],
                                                   background = background_s[,7:ncol(background_s)],
                                                   method = split_m),
                               silent = T)
                           
                           
                           
                         }
                         
                         
                         if(length(split_m) == 2){
                           
                           try(model_fold[[m]] <- 
                                 fit_plug_and_play(presence = presence_s[which(presence_data$fold!=fold),7:ncol(presence_s)],
                                                   background = background_s[,7:ncol(background_s)],
                                                   presence_method = split_m[1],
                                                   background_method = split_m[2]),
                               silent = T)
                           
                         }
                         
                         
                       } #m models loop end
                       
                       time_finish <- Sys.time()
                       
                       model_time <- time_finish - time_start
                       
                       #convert model time to seconds if needed
                       
                       if(units(model_time) != "secs"){ units(model_time) <- "secs" }
                       
                       if(units(model_time) != "secs"){stop("Model time units not seconds")}
                       
                       model_time <- as.numeric(model_time)
                       
                  ### end ensemble fitting
                       
                       
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
                       

            # Evaluate models
                       
                  training_predictions <- list()
                  testing_predictions <- list()
                  
            # Toss parts of model fold that weren't fitted
                  
                  
                  model_fold <-
                  model_fold %>%
                    discard(is.null)

                  
                for(m in 1:length(model_fold)){
                  
                  if("ratio" %in% names(model_fold[[m]])){
                    
                    training_predictions[[m]] <- project_density_ratio(dr_model = model_fold[[m]],
                                                                  data = training_data)
                    
                    testing_predictions[[m]] <- project_density_ratio(dr_model = model_fold[[m]],
                                                                 data = testing_data)

                  }else{
                    
                    training_predictions[[m]] <- project_plug_and_play(pnp_model = model_fold[[m]],
                                                                  data = training_data)
                    
                    testing_predictions[[m]] <- project_plug_and_play(pnp_model = model_fold[[m]],
                                                                 data = testing_data)
                    
                  }

                }#m loop for eval

            # end evaluate models
                  
                  training_predictions <- training_predictions %>%
                    as.data.frame() %>%
                    `colnames<-`(NULL)
              
                  testing_predictions <- testing_predictions %>%
                    as.data.frame() %>%
                    `colnames<-`(NULL)
                    
                  training_mean_vector = colMeans(training_predictions)
                  training_sd_vector = apply(X = training_predictions,MARGIN = 2,FUN = sd)
                  
                  std_testing_predictions <- (testing_predictions %>% as.matrix())/colSums(testing_predictions,na.rm = TRUE)
                  mean_testing_predictions <- rowMeans(std_testing_predictions,
                                                       na.rm=TRUE)
                  
                  # mean_testing_predictions <-
                  # testing_predictions %>%
                  # pbsdm:::rescale_w_objects(mean_vector = training_mean_vector,
                  #                           sd_vector = training_sd_vector) %>%
                  #   rowMeans(na.rm = TRUE)
                  
                  
                  std_training_predictions <- (training_predictions %>% as.matrix())/colSums(training_predictions,na.rm = TRUE)
                  mean_training_predictions <- rowMeans(std_training_predictions,
                                                       na.rm=TRUE)
                  
                  # 
                  # mean_training_predictions <-
                  #   training_predictions %>%
                  #   pbsdm:::rescale_w_objects(mean_vector = training_mean_vector,
                  #                             sd_vector = training_sd_vector) %>%
                  #   rowMeans(na.rm = TRUE)
                  # 
                  
                    

                       fold_training_suitability_v_occurrence <- data.frame(suitability = mean_training_predictions,
                                                                            occurrence = c(rep(1,length(which(presence_data$fold!=fold))),
                                                                                           rep(0,nrow(background_s))))
                       
                       fold_testing_suitability_v_occurrence <- data.frame(suitability = mean_testing_predictions,
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
                       
                       fold_threshold <- stats::quantile(x = fold_training_suitability_v_occurrence$suitability[which(fold_training_suitability_v_occurrence$occurrence==1)],
                                                    probs = quantile,
                                                    na.rm = T)
                       
                      # Code to make vote suitability
                       
                       for(r in 1:ncol(training_predictions)){
                         
                         threshold_r <- stats::quantile(x = training_predictions[,r][which(fold_training_suitability_v_occurrence$occurrence==1)],
                                                        probs = quantile,
                                                        na.rm = T)
                         if(r==1){
                           
                           testing_vote_suitability_out <- (testing_predictions[,r] >= threshold_r) %>%
                             as.numeric()
                           
                           training_vote_suitability_out <- (training_predictions[,r] >= threshold_r) %>%
                             as.numeric()
                           
                           
                         }else{
                           
                           testing_vote_suitability_out <- testing_vote_suitability_out +(testing_predictions[,r] >= threshold_r) %>%
                             as.numeric()
                           
                                                  
                           training_vote_suitability_out <- training_vote_suitability_out +(training_predictions[,r] >= threshold_r) %>%
                             as.numeric()
                           
                         }
                        
                        
                         
                          
                         
                       }
                       
                       # end code to make vote suitability
                       
                       vote_test_roc_obj <- tryCatch(pROC::roc(response = testing_vote_suitability_out,
                                                             predictor = fold_testing_suitability_v_occurrence$suitability,
                                                             level = c(0,1),
                                                             direction = "<"),
                                                   error = function(e){e})
                       
                       vote_training_roc_obj <- tryCatch(pROC::roc(response = training_vote_suitability_out,
                                                               predictor = fold_training_suitability_v_occurrence$suitability,
                                                               level = c(0,1),
                                                               direction = "<"),
                                                     error = function(e){e})
                       
                       
                       if(inherits(vote_test_roc_obj,"roc")){
                         
                         out$testing_vote_AUC[fold] <- vote_test_roc_obj$auc
                         
                         out$testing_vote_pAUC_specificity[fold] <- pROC::auc(roc = vote_test_roc_obj,
                                                                         partial.auc = c(.8, 1),
                                                                         partial.auc.correct = TRUE,
                                                                         partial.auc.focus = "specificity")[[1]]
                         
                         out$testing_vote_pAUC_sensitivity[fold] <- pROC::auc(roc = vote_test_roc_obj,
                                                                         partial.auc = c(.8, 1),
                                                                         partial.auc.correct = TRUE,
                                                                         partial.auc.focus = "sensitivity")[[1]]

                       }
                      
                       if(inherits(vote_training_roc_obj,"roc")){
                         
                         out$training_vote_AUC[fold] <- vote_training_roc_obj$auc
                         
                         out$training_vote_pAUC_specificity[fold] <- pROC::auc(roc = vote_training_roc_obj,
                                                                              partial.auc = c(.8, 1),
                                                                              partial.auc.correct = TRUE,
                                                                              partial.auc.focus = "specificity")[[1]]
                         
                         out$training_vote_pAUC_sensitivity[fold] <- pROC::auc(roc = vote_training_roc_obj,
                                                                              partial.auc = c(.8, 1),
                                                                              partial.auc.correct = TRUE,
                                                                              partial.auc.focus = "sensitivity")[[1]]
                         
                       }
                       
                      # vote threshold 
                       
                       vote_threshold <- stats::quantile(x = training_vote_suitability_out[which(fold_training_suitability_v_occurrence$occurrence == 1)],
                                                      probs = quantile,
                                                      na.rm = T)
                       

                       
                      # Anything greater than suitability threshold is considered a presence  
                       
                       TP <- length(which(fold_testing_suitability_v_occurrence$suitability >= fold_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 1))
                       
                       FN <- length(which(fold_testing_suitability_v_occurrence$suitability < fold_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 1))
                       
                       TN <- length(which(fold_testing_suitability_v_occurrence$suitability < fold_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 0))
                       
                       FP <- length(which(fold_testing_suitability_v_occurrence$suitability >= fold_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 0))
                       
                       # votes: Anything greater than suitability threshold is considered a presence  
                       
                       vote_TP <- length(which(testing_vote_suitability_out >= vote_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 1))
                       
                       vote_FN <- length(which(testing_vote_suitability_out < vote_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 1))
                       
                       vote_TN <- length(which(testing_vote_suitability_out < vote_threshold &
                                            fold_testing_suitability_v_occurrence$occurrence == 0))
                       
                       vote_FP <- length(which(testing_vote_suitability_out >= vote_threshold &
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
                       
                       vote_sensitivity <- vote_TP / (vote_TP + vote_FN)
                       vote_specificity <- vote_TN / (vote_FP + vote_TN)
                       vote_DOR <- (vote_TP*vote_TN)/(vote_FP*vote_FN)
                       #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
                       vote_prediction_accuracy <- (vote_TP+vote_TN)/(vote_TP+vote_TN+vote_FP+vote_FN)
                       vote_P_o <- (vote_TP+vote_TN)/(vote_TP+vote_TN+vote_FP+vote_FN)
                       vote_Ppres <- ((vote_TP+vote_FP)/(vote_TP+vote_TN+vote_FP+vote_FN))*((vote_TP+vote_FN)/(vote_TP+vote_TN+vote_FP+vote_FN))
                       vote_Pabs <- ((vote_FN+vote_TN)/(vote_TP+vote_TN+vote_FP+vote_FN))*((vote_FP+vote_TN)/(vote_TP+vote_TN+vote_FP+vote_FN))
                       vote_P_e <- vote_Ppres+vote_Pabs
                       vote_kappa <- (vote_P_o - vote_P_e)/(1-vote_P_e)
                       
                       
                       out$testing_DOR[fold] <- DOR
                       out$testing_prediction_accuracy[fold] <- prediction_accuracy
                       out$testing_sensitivity[fold] <- sensitivity
                       out$testing_specificity[fold] <- specificity
                       out$testing_kappa[fold] <- kappa
                       
                       out$testing_vote_DOR[fold] <- vote_DOR
                       out$testing_vote_prediction_accuracy[fold] <- vote_prediction_accuracy
                       out$testing_vote_sensitivity[fold] <- vote_sensitivity
                       out$testing_vote_specificity[fold] <- vote_specificity
                       out$testing_vote_kappa[fold] <- vote_kappa
                       
                       
                       
                       out$testing_correlation[fold] <- cor(fold_testing_suitability_v_occurrence$suitability,
                                                            fold_testing_suitability_v_occurrence$occurrence)
                       
                       out$testing_vote_correlation[fold] <- cor(testing_vote_suitability_out,
                                                                 fold_testing_suitability_v_occurrence$occurrence)
                       
                       
                       out$runtime[fold] <- model_time
                       
                       return(out %>%
                                rename(fold_temp = fold) %>%
                                filter(fold_temp == fold) %>%
                                rename(fold = fold_temp)
                       ) #renaming to prevent issues with having the same variable name in environment and dataframe
                       
                     }#end fold

      #Fit full model

      if(verbose){message("Starting full model fit")}
      
      time_start <- Sys.time()              
      
      model_full <- list()
      
      for(m in 1:length(model_vector)){
        
        model_m <- model_vector[m] 
        
        model_m <- gsub(pattern = " ",replacement = "",x = model_m)
        
        split_m <- strsplit(x = model_m,split = "/") %>% unlist()
        
        if(length(split_m) == 1){
          
          
          try(model_full[[m]] <- 
                fit_density_ratio(presence = presence_s[,7:ncol(presence_s)],
                                  background = background_s[,7:ncol(background_s)],
                                  method = split_m),
              silent = T)
          
        }
        
        
        if(length(split_m) == 2){
          
          try(model_full[[m]] <- 
                fit_plug_and_play(presence = presence_s[,7:ncol(presence_s)],
                                  background = background_s[,7:ncol(background_s)],
                                  presence_method = split_m[1],
                                  background_method = split_m[2]),
              silent = T)
          
        }
        
        
      }
      
      time_finish <- Sys.time()
      
      model_time_full <- time_finish - time_start
      
      # convert model time to seconds if needed
      
        if(units(model_time_full) != "secs"){ units(model_time_full) <- "secs" }
        
        if(units(model_time_full) != "secs"){stop("Model time units not seconds")}
        
        model_time_full <- as.numeric(model_time_full)
      
      # Full model projections
      
      if(verbose){message("Starting full model projections")}
        
      full_data <- rbind(presence_s[,7:ncol(presence_s)],
                         background_s[,7:ncol(background_s)])
      
      full_predictions <- list()

      # skip species if no models could be fit
      
      if(length(model_full)==0){
        message("Skipping model, was not fit")
        next
        }

      for(m in 1:length(model_full)){
        
        if(is.null(model_full[[m]])){next}
        
        if("ratio" %in% names(model_full[[m]])){
          
          full_predictions[[m]] <- project_density_ratio(dr_model = model_full[[m]],
                                                             data = full_data)
          
        }else{
          
          full_predictions[[m]] <- project_plug_and_play(pnp_model = model_full[[m]],
                                                             data = full_data)

        }
        
      } # m loop for eval

      
    # Make mean suitability score (re-scale relative suitability first)

      full_predictions <- full_predictions %>%
        discard(is.null) %>%
        as.data.frame() %>%
        `colnames<-`(NULL)
      
      std_full_predictions <- (full_predictions %>% as.matrix())/colSums(full_predictions,na.rm = TRUE)
      mean_full_predictions <- rowMeans(std_full_predictions,na.rm=TRUE)
      
# 
#       full_mean_vector <- colMeans(full_predictions)
#       
#       full_sd_vector <- apply(X = full_predictions,MARGIN = 2,FUN = sd)
#       
#       mean_full_predictions <-
#         full_predictions %>%
#         pbsdm:::rescale_w_objects(mean_vector = full_mean_vector,
#                                   sd_vector = full_sd_vector) %>%
#         rowMeans(na.rm = TRUE)
      
      
      
      
      full_suitability_v_occurrence <- data.frame(suitability = mean_full_predictions,
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
      
        full_threshold <- stats::quantile(x = full_suitability_v_occurrence$suitability[which(full_suitability_v_occurrence$occurrence==1)],
                                     probs = quantile,
                                     na.rm = T)
      
        #Anything greater than suitability threshold is considered a presence  
          
          TP <- length(which(full_suitability_v_occurrence$suitability >= full_threshold &
                               full_suitability_v_occurrence$occurrence == 1))
          
          FN <- length(which(full_suitability_v_occurrence$suitability < full_threshold &
                               full_suitability_v_occurrence$occurrence == 1))
          
          TN <- length(which(full_suitability_v_occurrence$suitability < full_threshold &
                               full_suitability_v_occurrence$occurrence == 0))
          
          FP <- length(which(full_suitability_v_occurrence$suitability >= full_threshold &
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
          pbsdm:::rescale_w_objects(data = pres_abs_data_s[,5:ncol(pres_abs_data_s)],
                            mean_vector = bg_means,
                            sd_vector = bg_sd)
      
  
        if(!all(pres_abs_data_s$siteid == data_i$pa$siteid[which(data_i$pa$spid == species_s)])){
          stop("Problem with data order in P/A data")
        }
      

      #Estimate suitabilities for PA data
            
        if(verbose){message("Starting estimation of PA ROR")}
        
        pa_predictions <- list()
        
        for(m in 1:length(model_full)){
          
          
          if(is.null(model_full[[m]])){next}
          
          
          if("ratio" %in% names(model_full[[m]])){
            
            pa_predictions[[m]] <- project_density_ratio(dr_model = model_full[[m]],
                                                           data = pres_abs_data_s[,5:ncol(pres_abs_data_s)])
            
          }else{
            
            pa_predictions[[m]] <- project_plug_and_play(pnp_model = model_full[[m]],
                                                           data = pres_abs_data_s[,5:ncol(pres_abs_data_s)])
            
          }
          
        }#m loop for eval
        
        
        # Make mean suitability score (re-scale relative suitability first)
            
            pa_predictions <- pa_predictions %>%
              discard(is.null) %>%
              as.data.frame() %>%
              `colnames<-`(NULL)
            
            std_pa_predictions <- (pa_predictions %>% as.matrix())/colSums(pa_predictions,na.rm = TRUE)
            mean_pa_predictions <- rowMeans(std_pa_predictions,na.rm=TRUE)
            
            # pa_mean_vector <- colMeans(pa_predictions)
            # 
            # pa_sd_vector <- apply(X = pa_predictions,MARGIN = 2,FUN = sd)
            # 
            # mean_pa_predictions <-
            #   pa_predictions %>%
            #   pbsdm:::rescale_w_objects(mean_vector = full_mean_vector,
            #                             sd_vector = full_sd_vector) %>%
            #   rowMeans(na.rm = TRUE)
            
            pa_suitability_v_occurrence <- 
              data.frame(suitability = mean_pa_predictions,
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
      
        full_threshold <- stats::quantile(x = full_suitability_v_occurrence$suitability[which(full_suitability_v_occurrence$occurrence==1)],
                                     probs = quantile,
                                     na.rm = T)
        
      #Anything greater than suitability threshold is considered a presence  
      
        TP <- length(which(pa_suitability_v_occurrence$suitability >= full_threshold &
                             pa_suitability_v_occurrence$occurrence == 1))
        
        FN <- length(which(pa_suitability_v_occurrence$suitability < full_threshold &
                             pa_suitability_v_occurrence$occurrence == 1))
        
        TN <- length(which(pa_suitability_v_occurrence$suitability < full_threshold &
                             pa_suitability_v_occurrence$occurrence == 0))
        
        FP <- length(which(pa_suitability_v_occurrence$suitability >= full_threshold &
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
        
      # code to evaluate model uncertainty via votes.
        # will consider two extremes: 1) any, 2) all
        
        if(verbose){message("Starting evaluation with any votes")}
        
        pa_all_vote_suitability_out <- NULL
        
        for(r in 1:ncol(full_predictions)){
          
          threshold_r <- stats::quantile(x = full_predictions[,r][which(full_suitability_v_occurrence$occurrence==1)],
                                         probs = quantile,
                                         na.rm = T)
          if(r==1){
            
            pa_all_vote_suitability_out <- (pa_predictions[,r] >= threshold_r) %>%
              as.numeric()
            
            
          }else{
            
            pa_all_vote_suitability_out <- pa_all_vote_suitability_out +(pa_predictions[,r] >= threshold_r) %>%
              as.numeric()
            
          }
        }
        
      # any models  
        
        
        TP <- length(which(pa_all_vote_suitability_out >= 1 &
                             pa_suitability_v_occurrence$occurrence == 1))
        
        FN <- length(which(pa_all_vote_suitability_out < 1 &
                             pa_suitability_v_occurrence$occurrence == 1))
        
        TN <- length(which(pa_all_vote_suitability_out < 1 &
                             pa_suitability_v_occurrence$occurrence == 0))
        
        FP <- length(which(pa_all_vote_suitability_out >= 1 &
                             pa_suitability_v_occurrence$occurrence == 0))
        
        out_full$pa_any_vote_sensitivity <- TP / (TP + FN)
        out_full$pa_any_vote_specificity <- TN / (FP + TN)
        #precision <- TP / (TP + FP)
        out_full$pa_any_vote_DOR <- (TP*TN)/(FP*FN)
        #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
        out_full$pa_any_vote_prediction_accuracy <- (TP+TN)/(TP+TN+FP+FN)
        
        P_o <- (TP+TN)/(TP+TN+FP+FN)
        Ppres <- ((TP+FP)/(TP+TN+FP+FN))*((TP+FN)/(TP+TN+FP+FN))
        Pabs <- ((FN+TN)/(TP+TN+FP+FN))*((FP+TN)/(TP+TN+FP+FN))
        P_e <- Ppres+Pabs
        out_full$pa_any_vote_kappa <- (P_o - P_e)/(1-P_e)
        
      # all models  
        
        if(verbose){message("Starting evaluation with all votes")}

        
        TP <- length(which(pa_all_vote_suitability_out >= ncol(full_predictions) &
                             pa_suitability_v_occurrence$occurrence == 1))
        
        FN <- length(which(pa_all_vote_suitability_out < ncol(full_predictions) &
                             pa_suitability_v_occurrence$occurrence == 1))
        
        TN <- length(which(pa_all_vote_suitability_out < ncol(full_predictions) &
                             pa_suitability_v_occurrence$occurrence == 0))
        
        FP <- length(which(pa_all_vote_suitability_out >= ncol(full_predictions) &
                             pa_suitability_v_occurrence$occurrence == 0))
        
        out_full$pa_all_vote_sensitivity <- TP / (TP + FN)
        out_full$pa_all_vote_specificity <- TN / (FP + TN)
        #precision <- TP / (TP + FP)
        out_full$pa_all_vote_DOR <- (TP*TN)/(FP*FN)
        #F1 <- 2*((precision * sensitivity)/(precision + sensitivity))
        out_full$pa_all_vote_prediction_accuracy <- (TP+TN)/(TP+TN+FP+FN)
        
        P_o <- (TP+TN)/(TP+TN+FP+FN)
        Ppres <- ((TP+FP)/(TP+TN+FP+FN))*((TP+FN)/(TP+TN+FP+FN))
        Pabs <- ((FN+TN)/(TP+TN+FP+FN))*((FP+TN)/(TP+TN+FP+FN))
        P_e <- Ppres+Pabs
        out_full$pa_all_vote_kappa <- (P_o - P_e)/(1-P_e)
        
        
        
      
      #Save output
      
        full_model_stats <- rbind(full_model_stats,
                                  data.frame(species = species_s, out_full))
        fold_model_stats <- rbind(fold_model_stats,
                                  data.frame(species = species_s, out))
        
        saveRDS(object = list(full_model_stats = full_model_stats,
                              fold_model_stats = fold_model_stats),
                file = temp_file)
        
        
      

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

