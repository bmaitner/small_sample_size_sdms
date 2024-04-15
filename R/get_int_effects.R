


get_int_effects <- function(model){

  if(inherits(x= model,what = "glmmTMB")){
    
    new_data <- data.frame(pres_method = unique(model$frame$pres_method)) %>%
      cross_join(data.frame(bg_method = unique(model$frame$bg_method))
      ) %>%
      mutate(n_presence = mean(model$frame$n_presence))
    
    response <- model$call$formula[[2]] %>% as.character()
    
  }
    
      
  if(inherits(model,"glm")){
    
    new_data <- data.frame(pres_method = unique(model$model$pres_method)) %>%
      cross_join(data.frame(bg_method = unique(model$model$bg_method))
      ) %>%
      mutate(n_presence = mean(model$model$n_presence))
    
    response <- model$formula[[2]] %>% as.character()

  }
  
  
  test <- predict(object = model,
                  new_data,
                  se.fit = TRUE,
                  type = "response")
  
  return( data.frame(response = response,
                     new_data,
                     test))  

}

