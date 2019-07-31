library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

# Script used in downstream_analysis/clinvar_data/preprocess_clinvar.Rmd workflow to parallelise the process of translating 
## DNA coordinates of ClinVar benign mutations to their protein coordinates. 
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

  clinvar.gr.benign <- clinvar.gr[clinvar.info$CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign") %>% unlist] # Filter for only pathogenic variants
  clinvar.gr.benign <- clinvar.gr.benign[seqnames(clinvar.gr.benign) %in%  c(1:22,"X","Y")] # filter for only autosomal and sex chromosomes
  # 126,086 variants are benign or likely benign 
  
  if (to > length(clinvar.gr.benign)){to == length(clinvar.gr.benign)}
  
  genome(clinvar.gr.benign) <- "GRCh37"
  benign.prot <- genomeToProtein(clinvar.gr.benign[from:to], EnsDb.Hsapiens.v75)
  benign.prot <- benign.prot %>% unlist
  write.csv(as.data.frame(benign.prot), sprintf("benign_prot_%s.csv", args[1]))
} 

main()