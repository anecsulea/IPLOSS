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

#################################################################################

for(sp in c("Chicken", "Duck")){

  load(paste(pathResults, sp, "/",annot,"/geneAnnotation.RData",sep=""))
  gene.ranges = geneAnnotation[["genes"]]

    for(sample in samples[[sp]]){
        path.res = paste(pathResults, sp, "/", annot, "/",sample,"/",sep="")

        load(paste(path.res,"matrix.after.clustering.RData",sep=""))

        system(paste("mkdir -p ",path.res,"peak_calling_clusters/",sep=""))

        clust = mtx.atac@meta.data$seurat_clusters

        for(cls in levels(clust)){

            print(cls)

            ## extract cells
            cells = rownames(mtx.atac@meta.data)[which(clust == cls)]
            writeLines(cells, con=paste(path.res,"peak_calling_clusters/cells_cluster",cls,".txt",sep=""))

            ## extract BAM

            system(paste("subset-bam --bam ",path.res, "outs/possorted_bam.bam --cell-barcodes ",path.res, "peak_calling_clusters/cells_cluster",cls,".txt --out-bam ", path.res, "peak_calling_clusters/cells_cluster",cls,".bam",sep=""))

            system(paste("macs2 callpeak -t ", path.res, "peak_calling_clusters/cells_cluster",cls,".bam -f BAMPE --outdir ",path.res, "peak_calling_clusters/ -n cluster",cls," -g 1.06e9 --keep-dup all --nomodel --shift 100 --extsize 200 --call-summits",sep=""))
      }

  }
}

#################################################################################

