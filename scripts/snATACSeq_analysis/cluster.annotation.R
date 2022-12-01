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

set.seed(1234)

#################################################################################

pathRNA="../../results/snRNASeq_analysis/"
pathResults="../../results/snATACSeq_analysis/"
annot="Ensembl"

samples=list("Chicken"=c("Gg_W_GT_male_33_scATAC"), "Duck"=c("Ap_W_GT_male_33_scATAC"))

#################################################################################

for(sp in c("Chicken", "Duck")){

    ## load gene annotation

    load(paste(pathResults, sp, "/",annot,"/geneAnnotation.RData",sep=""))

    gene.ranges = geneAnnotation[["genes"]]

    for(sample in samples[[sp]]){
        load(paste(pathResults, sp, "/", annot, "/",sample,"/matrix.afterQC.RData",sep=""))

        ## TF-IDF normalization

        DefaultAssay(mtx.filtered) <- 'bins'
        mtx.filtered = RunTFIDF(mtx.filtered)
        mtx.filtered = FindTopFeatures(mtx.filtered, min.cutoff = 'q99.9') ## bins with very high counts, we remove them afterwards
        top = VariableFeatures(mtx.filtered)

        ## not all bins are informative, so restrict to top 75%
        mtx.filtered = FindTopFeatures(mtx.filtered, min.cutoff = 'q25')
        print(length(VariableFeatures(mtx.filtered)))

        ## remove the ones that might be repetitive elements
        VariableFeatures(mtx.filtered)=VariableFeatures(mtx.filtered)[!(VariableFeatures(mtx.filtered) %in% top)]

        ## singular value decomposition
        mtx.filtered <- RunSVD(object = mtx.filtered,assay = 'bins',
                               reduction.key = 'LSI_',reduction.name = 'lsi')

        DepthCor(mtx.filtered) ## correlation between total read counts per cell and SVD components
        ## 1st component is highly correlated with total read count

        DimPlot(object = mtx.filtered, reduction = "lsi", dims = c(2,3),group.by = "orig.ident") + NoLegend()

        ## umap on components 2-30
        mtx.filtered <- RunUMAP(object = mtx.filtered, reduction = 'lsi', dims = 2:30)
        DimPlot(mtx.filtered, reduction = "umap")

        ## clustering
        mtx.filtered <- FindNeighbors(object = mtx.filtered, reduction = 'lsi', dims = 2:30)
        ## using default Louvain algorithm
        mtx.filtered <- FindClusters(object = mtx.filtered, verbose = FALSE, algorithm = 1, resolution = 0.8, random.seed = 1234)

        DimPlot(mtx.filtered, label = T, reduction = "umap")+ggtitle("res=0.8")

        #################################################################################
        ## link with scRNA-seq data

        ## load gene activity matrix computed with ArchR

        genematrix = read.table(paste0(pathResults, sp, "/", annot, "/",sample,"/ArchR_geneActivity.txt"), h=T, stringsAsFactors=F)
        colnames(genematrix) = unlist(lapply(colnames(genematrix), function(x) paste(unlist(strsplit(x, split="\\."))[-1], collapse="-")))

        ## select common cells
        filtered.cells = rownames(mtx.filtered@meta.data)
        common.cells = intersect(filtered.cells, colnames(genematrix))
        mtx.filtered = mtx.filtered[,which(colnames(mtx.filtered)%in%common.cells)]
        genematrix = genematrix[,rownames(mtx.filtered@meta.data)]

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
        mtx.filtered[['activity']] <- CreateAssayObject(data = genematrix)
    }
}

#################################################################################
