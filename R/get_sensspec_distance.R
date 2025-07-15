#code to get distance in sens-spec space
library(sf)
get_sensspec_distance <- function(combined_stats,
                                  self_comparison = TRUE){
  
  
  # get models
  
  models <- unique(combined_stats %>% pull(method))
  
  model_pairs <- combn(x = models,m = 2) %>% t()
  
  # If needed, self-comparisons should be added in here
  
  if(self_comparison){
    
    model_pairs <- rbind(model_pairs,cbind(models,models))
    
  }
  
  
  # iterate over pairs
  
    for(i in 1:nrow(model_pairs)){
      
      message("Model pair ",i," of ",nrow(model_pairs))
      
      pair_i <- model_pairs[i,]
      
      model_a <- pair_i[1]
      model_b <- pair_i[2]
      
      
      # pull the relevant data for all species/sites and both models
      
      data_a <- combined_stats %>%
        filter(method == model_a) %>%
        collect() %>%
        select(species,
               sens_a = pa_sensitivity,
               spec_a = pa_specificity)#not sure if its better to use training or testing sens/spec

      data_b <- combined_stats %>%
        filter(method == model_b) %>%
        collect() %>%
        select(species,,
               sens_b = pa_sensitivity,
               spec_b = pa_specificity) #not sure if its better to use training or testing sens/spec

      
      #join data
      
      combined_data <- full_join(x = data_a,
                y = data_b,
                by = "species")
      

      rm(data_a,data_b)
      
      # summarize into the needed stats  
      
        combined_data <- 
        combined_data %>%
          mutate(sens_spec_distance = sqrt(((sens_a - sens_b)^2) +
                                             ((spec_a - spec_b)^2)
                                           )
                 ) %>%
          select(species,sens_spec_distance) %>%
          mutate(model_a = model_a,
                 model_b = model_b)
        
        
        if(i == 1){
          
          out <- combined_data
          
        }else{
          
          out <- bind_rows(out,combined_data)
          
        }
        
      
      
    }#i loop
  
  
  return(out)    
  
  
}