# can only use AUC, because thresholding not done in Valavi paper.

library(tidyverse)
library(disdat)
get_valavi_stats <- function(models_prediction_folder = "data/manual_downloads/Valavi/Models_prediction/"){
  
  # make output file
  output <- NULL
  
  # list folders
  
  model_folders <- list.dirs(models_prediction_folder,full.names = TRUE) %>%
    setdiff(models_prediction_folder)

  
  for(i in 1:length(model_folders)){
    
    dir_i <- model_folders[i]
    files_i <- list.files(dir_i,recursive = TRUE,full.names = TRUE)
    
    for(j in 1:length(files_i)){
      
      # get data
      
      data_j <- read.csv(files_i[j])

      #join with PA data
      
      data_j <- disdat::disData(data_j$region)$pa %>%
        dplyr::filter(spid == unique(data_j$spid)) %>%
        select(spid,siteid,pa)%>%
        right_join(data_j )
      
      # calc AUC
      
          training_roc_obj <- tryCatch(pROC::roc(response = data_j$pa,
                                                 predictor = data_j$prediction,
                                                 level = c(0,1),
                                                 direction = "<"),
                                       error = function(e){e})
          
          # only attempt to use ROC if it could be calculated properly
          
            if(inherits(training_roc_obj,"roc")){
              
              
              AUC_j <- training_roc_obj$auc
      
            }else{
            
              AUC_j <- NA
              
            }
          
          
          data_ij <-
          data_j %>%
            select(spid,region,model)%>%
            unique()%>%
            mutate(pa_AUC = AUC_j %>% as.numeric(),
                   pa_correlation = cor(x = data_j$prediction,
                                         y = data_j$pa,
                                         method = "pearson")
                   )
          
          output <- bind_rows(output,data_ij)
          
          
    }#j loop
    
  }# end i loop

return(output)  
  
}# fx