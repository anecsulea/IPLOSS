#################################################################################

## adapted from original script by Menghan Wang

#################################################################################

library(Signac)
library(Seurat)
library(ArchR)
library(dplyr)
library(GenomeInfoDb)
library(ensembldb)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(GenomicRanges)
library(future)
library(pheatmap)

set.seed(1234)

#################################################################################

pathResults="../../results/snATACSeq_analysis/"
annot="Ensembl"

samples=list("Chicken"=c("Gg_W_GT_male_33_scATAC"), "Duck"=c("Ap_W_GT_male_33_scATAC"))

samplesRNA=c("Gg_W_GT_male_33_scRNA", "Ap_W_GT_male_33_scRNA")
names(samplesRNA)=c("Gg_W_GT_male_33_scATAC", "Ap_W_GT_male_33_scATAC")

#################################################################################

for(sp in c("Chicken", "Duck")){
   
  load(paste(pathResults, sp, "/",annot,"/geneAnnotation.RData",sep=""))
  gene.ranges = geneAnnotation[["genes"]]
  
  for(sample in samples[[sp]]){
    load(paste(pathResults, sp, "/", annot, "/",sample,"/matrix.after.clustering.RData",sep=""))
  }
}

#################################################################################
    
