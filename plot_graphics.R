# A script to generate clustering plots and quality score plot from results of spectrus
library(dplyr)
library(ggplot2)

main <- function(){
  # Accepting Command line Argument
  args <- commandArgs(trailingOnly = TRUE)
  dir.name <- args[1]
  output.clustering <- args[2]
  output.quality <- args[3]
 
  # Plotting clustering of a protein sequence
  ## read every cluster labels in a directory 
  name.array <- list.files(dir.name, pattern = 'final_labels_kmed-*')
  df <- data.frame()
  for (i in 1:length(name.array)){
    df.single <- read.csv(paste(dir.name, name.array[i], sep = "/"), header = FALSE, sep = " ")
    df.single$V3 <- gsub(".dat", "", gsub("final_labels_kmed-", "", name.array[i])) %>% as.numeric()
    df.single$V2 <- df.single$V2 %>% as.factor
    colnames(df.single) <- c("residue_position", "membership", "q")
    df = rbind(df, df.single) 
  }
  df = df[order(df$q,df$residue_position),]
  df = df[df$q <= 10,] #only plot the top 10 
  
  ## define coloring scheme for 10 clusterings 
  #### TODO: make this variable dynamic to the number of clusters detected #### 
  colors = c("#4e79a7", "#59a14f", "#9c755f", "#f28e2b", "#edc948", "#e15759", "#b07aa1",
             "#76b6b2", "#ff9da7", 	"#000000", "#bab0ac")
  
  ## use scatter plot with square dots
  df1 = df; df1$q <- df$q + .05
  df2 = df; df2$q <- df$q + .1
  df3 = df; df3$q <- df$q + .15
  df4 = df; df4$q <- df$q + .2
  
  df.combined = rbind(df, df1, df2, df3, df4)
  
  colnames(df.combined) <- c("residue_position", "membership", "q")
  
  
  options(bitmapType='cairo') 
  png(output.clustering)
  print(ggplot(df.combined, aes(x=residue_position, y=q, color=membership)) + 
    geom_point(size=1.3, shape=15) + 
    scale_color_manual(values=colors) + labs(x="residue position", y="number of clusters"))
  dev.off()
  

  # Plot quality score against number of clusters
  system(sprintf("sed -i 's/^[ \t]*//' %squality_score.dat", dir.name))
  df.qualscore <- read.csv(paste(dir.name, "quality_score.dat", sep = "/"), header = FALSE, strip.white=TRUE, sep = " ")
  colnames(df.qualscore) <- c("q", "median_quality_score", "mean_quality_score", "prefactor1", "prefactor2")
  df.qualscore = df.qualscore[df.qualscore$q <= 10,]
  
  png(output.quality)
  print(ggplot(df.qualscore, aes(x=q, y=median_quality_score)) + geom_line() + geom_point() + 
    labs(x="number of clusters", y="quality score"))
  dev.off()
}

main()
