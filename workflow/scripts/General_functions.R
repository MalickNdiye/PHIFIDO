parse_sample_name<- function(x){
  name_list<- unlist(strsplit(x, split = "-"))
  
  strain=name_list[1]
  cocktail=name_list[2]
  batch=name_list[3]
  
  out<- c("strain"=strain, "cocktail"=cocktail, "batch"=batch)
  
  return(out)
}

get_mtadata<- function(df, toget, corre, corre_values){
  # This function returns a given "toget" data from a table
  toget <- as.character(toget)
  corre <- as.character(corre)
  
  # Check columns exist
  if (!toget %in% names(df)) {
    stop(sprintf("Column '%s' not found in dataframe", toget))
  }
  if (!corre %in% names(df)) {
    stop(sprintf("Column '%s' not found in dataframe", corre))
  }
  
  # Create a named lookup vector
  lookup <- setNames(df[[toget]], df[[corre]])
  
  # Match corre_values to lookup names
  result <- lookup[as.character(corre_values)]
  
  # Always return a plain vector
  as.vector(result)
  
}

Bifido_strain_order<- factor(c("ESL0820", "ESL1077", "ESL0819", "ESL0200", "ESL0824", "ESL0198", "ESL0827", "ESL0199", "ESL0822", "ESL0170", "ESL1069", "ESL1060"),
                             levels= c("ESL0820", "ESL1077", "ESL0819", "ESL0200", "ESL0824", "ESL0198", "ESL0827", "ESL0199", "ESL0822", "ESL0170", "ESL1069", "ESL1060"))
