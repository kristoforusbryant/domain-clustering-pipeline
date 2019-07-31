library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

# Script used in downstream_analysis/clinvar_data/preprocess_clinvar.Rmd workflow to parallelise the process of translating 
## DNA coordinates of ClinVar pathogenic mutations to their protein coordinates. 
main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  gap <- 10000 
  from <- (as.numeric(args[1]) - 1) * gap  + 1 
  to <- from + gap - 1 
  
  clinvar.vcf <- readVcf("clinvar_20190715.vcf", genome = "hg19")
  clinvar.vcf %>% header 
  clinvar.gr <- rowRanges(clinvar.vcf)
  clinvar.info <- clinvar.vcf %>% info()
  # total 453,467 annotated variants
  
  clinvar.gr.patho <- clinvar.gr[clinvar.info$CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic") %>% unlist] # Filter for only pathogenic variants
  clinvar.gr.patho <- clinvar.gr.patho[seqnames(clinvar.gr.patho) %in%  c(1:22,"X","Y")] # filter for only autosomal and sex chromosomes
  # 92,822 variants are pathologic or likely pathologic
  
  if (to > length(clinvar.gr.patho)){to == length(clinvar.gr.patho)}
  
  genome(clinvar.gr.patho) <- "GRCh37"
  pathogenic.prot <- genomeToProtein(clinvar.gr.patho[from:to], EnsDb.Hsapiens.v75) # Run in parallel on server
  pathogenic.prot <- pathogenic.prot %>% unlist
  write.csv(as.data.frame(pathogenic.prot), sprintf("pathogenic_prot_%s.csv", args[1]))
} 

main()