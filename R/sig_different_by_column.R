#assumes first column is grouping variable, others are to be compared
sig_different_by_column <- function(output_df,filter_by = "pa_AUC"){
  
  #make output dfs
  
  output_df %>%
    select(1)%>%
    unique()->model_types

  df<-matrix(nrow = nrow(model_types),
         ncol = ncol(output_df))%>%as.data.frame()
  
  colnames(df) <- colnames(output_df)
  df[,1] <- model_types[,1]
  
  if(!is.null(filter_by)){
    
    df %>%
      relocate(1,all_of(filter_by))->df
    
  }
  
  W_table <- df
  p_val_table <- df  

  #loop over rows
  
  for(i in 2:ncol(df)){
    
    #variable 
    grouping_var <- colnames(df)[1]  
    var_i <- colnames(df)[i]
    
    #get data subset
    
    output_df_i <-
      output_df %>%
      select(all_of(c(grouping_var,var_i)))%>%
      rename(model= 1,
             metric = 2)
    
    # get best performing
    
    if(!is.null(filter_by)){
      if(i>2){
        
        #first filter out anything with a fitlering value worse than the best
        
        equivalent_to_best <-
          p_val_table%>%
          select(1,all_of(filter_by))%>%
          rename(filter_by = 2) %>%
          filter(filter_by > 0.05) %>%
          pull(model)
        
        
        best_model_i <-
          output_df_i %>%
          group_by(model)%>%
          filter(model %in% equivalent_to_best)%>%
          summarise(mean = mean(metric,na.rm = TRUE)) %>%
          arrange(-mean) %>%
          slice_head(n = 1) %>%
          pull(model)  
        
        
      }else{
        
        best_model_i <-
          output_df_i %>%
          group_by(model)%>%
          summarise(mean = mean(metric,na.rm = TRUE)) %>%
          arrange(-mean) %>%
          slice_head(n = 1) %>%
          pull(model)
        
        
      }
    }else{
      
      best_model_i <-
        output_df_i %>%
        group_by(model)%>%
        summarise(mean = mean(metric,na.rm = TRUE)) %>%
        arrange(-mean) %>%
        slice_head(n = 1) %>%
        pull(model)
      
      
    }
    
    
    
    
    # get distribution of best model
    
    data_y <- output_df_i %>%
      filter(model==best_model_i) %>%
      pull(metric)
    
    
    # loop over and compare with all models
    
    for(j in 1:nrow(df)){
      
      data_x <- output_df_i %>%
        filter(model == df$model[j]) %>%
        pull(metric)
    
      ti <-tryCatch(wilcox.test(x = data_x,
                                y = data_y),
                    error=function(e){e})
      
      
      if(inherits(ti,"error")){
        ti$p.value <- NA
        ti$statistic <- NA
      }
      
      p_val_table[j,i]<- round(x = ti$p.value,digits =  4) 
      W_table[j,i]<- ti$statistic
      
      
    }#j
    
    

    
  }#i
  
return(list(p_val_table = p_val_table,
       W_table = W_table))
  
  
    
}
