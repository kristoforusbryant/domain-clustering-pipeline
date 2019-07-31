# For final data containing all variants, summarise the data into results dataframe
# containing ENSG, cluster_membership, n_of_mutation_in_that_domain, 
# n_of_residue_in_that_domain
library(dplyr) 

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  gap = 1000
  from = (as.numeric(args[1]) - 1) * gap + 1  
  to = from + gap -1 
  
  print("Starting...")
  print(from) 
  print(to)
   
  df <- read.csv("final_all_variants.csv") # 2,862,793
  df$AF <- as.character(df$AF) %>% as.numeric()
  df[is.na(df$AF),]$AF <- 0. #1,493,667 rows with no mutations 
  
  result.df <- df[,c(1,3)] %>% unique() 
  result.df <- result.df[order(result.df$ENSG, result.df$cluster_membership),]
  if (to > nrow(result.df)){to = nrow(result.df)}
  result.df.reduced <- result.df[from:to,] 
  result.df.reduced$cluster_size <- sapply(from:to, function(i){ df[df$ENSG == result.df[i,1] & 
                                                                        df$cluster_membership == result.df[i,2],][c(1,2)] %>%  
									unique %>% nrow})
  
  result.df.reduced$n.mut <- sapply(from:to, function(i){ df[df$ENSG == result.df[i,1] &  
                                                                 df$cluster_membership == result.df[i,2] & 
                                                                 df$AF > 0,] %>% nrow})
  
  write.csv(result.df.reduced, sprintf("final_reduced_%s.csv", args[1]), row.names = FALSE)
  print("DONE") 
}

main()
