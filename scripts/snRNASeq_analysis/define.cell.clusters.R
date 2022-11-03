## original script by Menghan Wang, University of Basel

###########################################################################

pathResults="../../results/snRNASeq_analysis/"
pathAnnot="../../data/ensembl_annotations/"
pathCellCycle="../../data/single_cell_data/cell_cycle_genes/"

splist=c("Chicken", "Duck")
annot=c("Ensembl", "EnsemblStringTie")
samples=list("Chicken"=c("Gg_W_GT_male_33_scRNA"), "Duck"=c("Ap_W_GT_female_33_scRNA",  "Ap_W_GT_male_33_scRNA"))
sample.names=c("gg_male", "ap_female", "ap_male")
names(sample.names)=c("Gg_W_GT_male_33_scRNA", "Ap_W_GT_female_33_scRNA",  "Ap_W_GT_male_33_scRNA")

###########################################################################

## load cell cycle genes

load(paste(pathCellCycle, "chicken_cycle_markers.Rda",sep=""))
load(paste(pathCellCycle, "duck_cycle_markers.RDa",sep=""))

cc.pairs=list("Chicken"=gg.pairs, "Duck"=ap.pairs)

###########################################################################

## load libraries

library(Seurat)
library(scran)
library(BiocParallel)

set.seed(42)

###########################################################################

for(sp in splist){
    gene.names=read.table(paste(pathAnnot, sp, "/GeneNames_Ensembl103.txt",sep=""),h=T, stringsAsFactors=F,sep="\t")
    rownames(gene.names)=gene.names[,1]

    for(sample in samples[[sp]]){
        for(annot in c("Ensembl", "EnsemblStringTie")){

            print(paste(sp, sample, annot))

            short.name=sample.names[sample]

            ## save all plots here
            pdf(file=paste(pathResults, sp, "/",annot, "/",sample, "/Clustering.pdf",sep=""))

            ## load previously created Seurat object

            load(paste(pathResults, sp, "/",annot, "/",sample, "/seurat.object.afterQC.RData",sep=""))

            ## rename cells to add short sample name

            newnames=paste0(short.name,"_",rownames(filtered.seurat@meta.data))
            filtered.seurat = RenameCells(filtered.seurat ,new.names = newnames)

            all.genes = rownames(filtered.seurat@assays$RNA@counts)

            ## log-normalize count data
            filtered.seurat = NormalizeData(filtered.seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)

            ## center the data (no scaling)
            filtered.seurat = ScaleData(filtered.seurat, features = all.genes, verbose = F, do.scale = F)

            ## take into account cell cycle
            ## use only S and G2/M phases
            pairs.sm = cc.pairs[[sp]][-1]

            ## Assign scores using the scaled data
            assign.se = cyclone(filtered.seurat@assays$RNA@scale.data, pairs.sm, gene.names=rownames(filtered.seurat), BPPARAM=bpparam(), verbose=T)

            ## The difference between the two phases
            rownames(assign.se$scores) = colnames(filtered.seurat)
            assign.se$scores$CC.Difference = (assign.se$scores$S - assign.se$scores$G2M)

            #Add the cc difference to the Seurat object
            filtered.seurat = AddMetaData(filtered.seurat, metadata = assign.se$scores)

            ## change type for percent.mt info
            filtered.seurat[["percent.mt"]] =  filtered.seurat$percent.mt

            ## Transform the data, removing sources of variation
            filtered.seurat = SCTransform(filtered.seurat, vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"),
                                          verbose = FALSE, return.only.var.genes = F, variable.features.n = NULL)

            detach("package:scran", unload=TRUE)

            ## extract variable features
            var.features=HVFInfo(filtered.seurat, assay="SCT")

            # only those with a variability (here residual variance)  > median + MAD
            var.features = rownames(var.features)[which(var.features[,3] > ( median(var.features[,3]) + mad(var.features[,3])))]

            ## ? should we remove W chromosome here ?

            ## Scale the original data as well, to show differences - is this necessary?
            ## filtered.seurat = ScaleData(filtered.seurat, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = F)

            ## run PCA on SCT data

            filtered.seurat=RunPCA(filtered.seurat, assay = "SCT", verbose = F, features = var.features)
            PCAPlot(filtered.seurat, group.by="orig.ident")

            ## eigenvalues
            barplot(filtered.seurat@reductions$pca@stdev^2)

            ## we can keep 30 dimensions, probably more than needed
            filtered.seurat=RunUMAP(filtered.seurat, reduction = "pca", dims = 1:30, verbose = F)

            ## umap plot
            DimPlot(filtered.seurat, reduction="umap")


            ## do the actual clustering
            filtered.seurat=FindNeighbors(filtered.seurat, dims=1:30, verbose = F)

            ## we go for lower resolution, too much detail in the mesenchyme cluster
            filtered.seurat=FindClusters(filtered.seurat, resolution = 0.2, algorithm = 4, verbose = F, random.seed = 42)

            clusters=attr(filtered.seurat, "meta.data")$seurat_clusters

            ## find markers

            markers=list()

            for(i in levels(clusters)){
                this.markers = FindMarkers(filtered.seurat, ident.1 = i, min.pct = 0.25)
                this.markers = this.markers[which(this.markers$avg_log2FC > 0.5 & this.markers$p_val_adj<0.01),]

                this.markers$GeneName=gene.names[rownames(this.markers),2]

                markers[[i]]=this.markers
            }

            save(list=c("filtered.seurat", "markers"), file=paste(pathResults, sp, "/",annot, "/",sample, "/seurat.object.after.clustering.RData",sep=""))

            ## close pdf
            dev.off()
        }
  }
}

###########################################################################
