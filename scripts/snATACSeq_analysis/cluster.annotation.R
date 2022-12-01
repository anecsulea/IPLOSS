#################################################################################

## adapted from original script by Menghan Wang

#################################################################################

library(Signac)
library(Seurat)
library(dplyr)
library(GenomeInfoDb)
library(ensembldb)
library(ggplot2)
library(patchwork)
library(scales)

set.seed(1234)

#################################################################################

pathRNA="../../results/snRNASeq_analysis/"
pathResults="../../results/snATACSeq_analysis/"
annot="Ensembl"

samples=list("Chicken"=c("Gg_W_GT_male_33_scATAC"), "Duck"=c("Ap_W_GT_male_33_scATAC"))

samplesRNA=c("Gg_W_GT_male_33_scRNA", "Ap_W_GT_male_33_scRNA")
names(samplesRNA)=c("Gg_W_GT_male_33_scATAC", "Ap_W_GT_male_33_scATAC")

#################################################################################

for(sp in c("Chicken", "Duck")){

  ## load gene annotation
  
  load(paste(pathResults, sp, "/",annot,"/geneAnnotation.RData",sep=""))
  
  gene.ranges = geneAnnotation[["genes"]]
  
  for(sample in samples[[sp]]){
    load(paste(pathResults, sp, "/", annot, "/",sample,"/matrix.afterQC.RData",sep=""))
    
    mtx.atac = mtx.filtered
    
    ## TF-IDF normalization
    
    DefaultAssay(mtx.atac) <- 'bins'
    mtx.atac = RunTFIDF(mtx.atac)
    mtx.atac = FindTopFeatures(mtx.atac, min.cutoff = 'q99.9') ## bins with very high counts, we remove them afterwards
    top = VariableFeatures(mtx.atac)
    
    ## not all bins are informative, so restrict to top 75%
    mtx.atac = FindTopFeatures(mtx.atac, min.cutoff = 'q25')
    print(length(VariableFeatures(mtx.atac)))
    
    ## remove the ones that might be repetitive elements
    VariableFeatures(mtx.atac)=VariableFeatures(mtx.atac)[!(VariableFeatures(mtx.atac) %in% top)]
    
    ## singular value decomposition
    mtx.atac <- RunSVD(object = mtx.atac,assay = 'bins',
                       reduction.key = 'LSI_',reduction.name = 'lsi')
    
    DepthCor(mtx.atac) ## correlation between total read counts per cell and SVD components
    ## 1st component is highly correlated with total read count
    
    DimPlot(object = mtx.atac, reduction = "lsi", dims = c(2,3),group.by = "orig.ident") + NoLegend()
    
    ## umap on components 2-30
    mtx.atac <- RunUMAP(object = mtx.atac, reduction = 'lsi', dims = 2:30)
    DimPlot(mtx.atac, reduction = "umap")
    
    ## clustering
    mtx.atac <- FindNeighbors(object = mtx.atac, reduction = 'lsi', dims = 2:30)
    ## using default Louvain algorithm
    mtx.atac <- FindClusters(object = mtx.atac, verbose = FALSE, algorithm = 1, resolution = 0.8, random.seed = 1234)
    
    DimPlot(mtx.atac, label = T, reduction = "umap")+ggtitle("res=0.8")
    
   #################################################################################
    ## link with scRNA-seq data
    
    ## load gene activity matrix computed with ArchR
    
    genematrix = read.table(paste0(pathResults, sp, "/", annot, "/",sample,"/ArchR_geneActivity.txt"), h=T, stringsAsFactors=F)
    colnames(genematrix) = unlist(lapply(colnames(genematrix), function(x) paste(unlist(strsplit(x, split="\\."))[-1], collapse="-")))
    
    ## select common cells
    atac.cells = rownames(mtx.atac@meta.data)
    common.cells = intersect(atac.cells, colnames(genematrix))
    mtx.atac = mtx.atac[,which(colnames(mtx.atac)%in%common.cells)]
    genematrix = genematrix[,rownames(mtx.atac@meta.data)]
    
    ## remove genes with ambiguous names
    geneinfo = data.frame(gene_id=gene.ranges$gene_id,gene_name=gene.ranges$gene_name)
    geneinfo = geneinfo[which(geneinfo$gene_name %in% rownames(genematrix)),]
    dupli = geneinfo$gene_name[which(duplicated(geneinfo$gene_name))]
    geneinfo = geneinfo[which(!geneinfo$gene_name%in%dupli),]
    rownames(geneinfo)=geneinfo$gene_name
    
    genematrix = genematrix[which(rownames(genematrix)%in%geneinfo$gene_name),]
    rownames(genematrix)=geneinfo[rownames(genematrix), "gene_id"] ## Ensembl gene ids
    
    ## log transform
    genematrix=as.matrix(genematrix)
    genematrix = as(log(genematrix+1), "dgCMatrix")
    
    ## add to object
    mtx.atac[['activity']] <- CreateAssayObject(data = genematrix)
    
    ## load scRNA-seq data
    
    sampleRNA = samplesRNA[sample]
    
    load(paste(pathRNA, sp, "/", annot, "/", sampleRNA, "/seurat.object.after.clustering.RData",sep=""))
    
    mtx.rna = atac.seurat
    DefaultAssay(mtx.rna) <- 'RNA'
    
    ## get variable features
    mtx.rna = FindVariableFeatures(mtx.rna)
    hvgs = VariableFeatures(mtx.rna)
    hvgs = intersect(hvgs, rownames(mtx.atac[["activity"]]))
    
    print(table(Idents(mtx.rna)))
    
    mtx.rna@meta.data$annotation = Idents(mtx.rna)
    
    ## as the cell type annotation is not done yet, use the cluster id to run following:
    mtx.rna@meta.data$annotation = paste0("ctype",Idents(mtx.rna))
    
    ## label transferring
    
    transfer.anchors = FindTransferAnchors(
      reference = mtx.rna,
      query = mtx.atac,
      dims = 1:30,
      reference.assay = "RNA", query.assay = "activity",
      max.features = 200, k.filter = 200,
      reduction = 'cca', features = hvgs
      )
    
    predicted.labels = TransferData(
      anchorset = transfer.anchors,
      refdata = Idents(mtx.rna),
      weight.reduction = mtx.atac[["lsi"]],
      dims = 1:30 ## same as for transfer anchors. This doesn't change results significantly
      )
    
    mtx.atac = AddMetaData(mtx.atac, metadata = predicted.labels)
    thres = 0.5
    mtx.atac@meta.data[mtx.atac@meta.data$prediction.score.max < thres,"predicted.id"]="NA"
    DimPlot(mtx.atac, reduction = "umap",group.by = "predicted.id")

    save(mtx.atac, file=paste(pathResults, sp, "/", annot, "/",sample,"/matrix.after.clustering.RData",sep=""))
  }
}

#################################################################################
