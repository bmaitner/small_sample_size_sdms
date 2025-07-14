#code to get the overlap in predicted locations
library(arrow)
library(tidyverse)

get_prediction_overlap <- function(dataset_folder = "outputs/model_predictions/",
                                   self_comparison = TRUE){
  
  
  dataset_info <- open_dataset(dataset_folder)  
  
  # dataset_info
    # siteid: string
    # group: string
    # species: string
    # quantile: double
    # predicted_occurrence: double
    # model: string
  
  # Get models
  
    models <- dataset_info %>%
      select(model) %>%
      unique() %>%
      collect() %>%
      pull(model)
    
  # Get species
    
    spp <- dataset_info %>%
      select(species) %>%
      unique() %>%
      collect() %>%
      pull(species)
  
  
  # Get pairwise stuff to iterate over
    
      model_pairs <- combn(x = models,m = 2) %>% t()
      
  # If needed, self-comparisons should be added in here
      
      if(self_comparison){
        
        model_pairs <- rbind(model_pairs,cbind(models,models))

      }

  # iterate over pairs
      
      for(i in 1:nrow(model_pairs)){
        
        pair_i <- model_pairs[i,]
        
        model_a <- pair_i[1]
        model_b <- pair_i[2]
        
        # pull the relevant data for all species/sites and both models
        
          data_a <- dataset_info %>%
                      filter(model == model_a) %>%
            collect() %>%
            mutate(model = "a")
  
          data_b <- dataset_info %>%
            filter(model == model_b) %>%
            collect() %>%
            mutate(model = "b")
        
        #join data
          
          
          combined_data <-
          data_a %>%
            bind_rows(data_b) %>%
            group_by(siteid,group,species,quantile)%>%
            pivot_wider(names_from = model,
                        values_from = predicted_occurrence)
          
          rm(data_a,data_b)
    
        # summarize into the needed stats  
          # model agreement (# agreed/ # total)
          # Jaccard (intersection/union) #gives agreement based on presences only
          
          combined_data <-
          combined_data %>%
            group_by(species,quantile) %>%
            mutate(agreement = as.numeric(a==b),
                   agreed_present = as.numeric(a== 1 & b == 1),
                   one_or_more_present = as.numeric(a==1 | b ==1)) %>%
            summarise(model_agreement = sum(agreement)/n(),
                      Jaccard = sum(agreed_present)/sum(one_or_more_present)) %>%
            mutate(model_a = model_a,
                   model_b = model_b)
            
        if(i == 1){
          
          out <- combined_data
          
        }else{
          
          out <- bind_rows(out,combined_data)
          
        }
        

      }# i loop
  

  return(out)    
  
}
