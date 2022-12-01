#################################################################################

## adapted from original script by Menghan Wang

#################################################################################

library(ArchR)
library(Seurat)
library(dplyr)
library(BSgenome)
library(GenomeInfoDb)
library(ensembldb)
library(ggplot2)
library(patchwork)

set.seed(1234)

addArchRThreads(threads = 4)

#################################################################################

pathResults="../../results/snATACSeq_analysis/"
annot="Ensembl"

samples=list("Chicken"=c("Gg_W_GT_male_33_scATAC"), "Duck"=c("Ap_W_GT_male_33_scATAC"))

#################################################################################

for(sp in c("Chicken", "Duck")){

    path.annot = paste(pathResults, sp, "/", annot, "/", sep="")

    load(paste0(path.annot,"geneAnnotation.RData"))
    load(paste0(path.annot,"genomeAnnotation.RData"))


    for(sample in samples[[sp]]){

        ## Create Arrow Files & ArchRProject

        ArrowFiles = createArrowFiles(
            inputFiles = paste0(pathResults, sp, "/", annot, "/",sample,"/outs/fragments.tsv.gz"),
            sampleNames = sample,
            geneAnnotation = geneAnnotation,
            genomeAnnotation = genomeAnnotation,
            minFrags = 100,
            minTSS = 1,
            addTileMat = TRUE,
            addGeneScoreMat = TRUE,
            force = TRUE
        )

        proj = ArchRProject(
            ArrowFiles = ArrowFiles,
            outputDirectory = paste0(pathResults, sp, "/", annot, "/",sample),
            geneAnnotation = geneAnnotation,
            genomeAnnotation = genomeAnnotation,
            copyArrows = FALSE
        )

        ## get gene score matrix

        tmp = getMatrixFromProject(proj,useMatrix = "GeneScoreMatrix")
        genematrix = tmp@assays@data[[1]] ## dgcMatrix
        rownames(genematrix)= tmp@elementMetadata$name ## genenames of dgcMatrix

        ## save data

        write.table(genematrix,paste0(pathResults, sp, "/", annot, "/",sample,"/ArchR_geneActivity.txt"), col.names = TRUE, row.names = TRUE, sep = "\t")

        saveArchRProject(ArchRProj = proj, outputDirectory=paste0(pathResults, sp, "/", annot, "/",sample,"/"))

        ## doublet assignment

        doubScores <- addDoubletScores(
            input = paste0(pathResults, sp, "/", annot, "/",sample,"/ArrowFiles/",sample,".arrow"),
            k = 10,
            knnMethod = "UMAP",
            LSIMethod = 1
        )

        ## doublet info was saved to disk, need to re-read project

        proj <- ArchRProject(
            ArrowFiles = paste0(pathResults, sp, "/", annot, "/",sample,"/ArrowFiles/",sample,".arrow"),
            outputDirectory = paste0(pathResults, sp, "/", annot, "/",sample),
            geneAnnotation = geneAnnotation,
            genomeAnnotation = genomeAnnotation,
            copyArrows = FALSE
        )

        projfiltered = filterDoublets(proj)
        allcells = rownames(getCellColData(proj))
        doublets = allcells[!(allcells %in% rownames(getCellColData(projfiltered)))]

        write.table(doublets,paste0(pathResults, sp, "/", annot, "/",sample,"/ArchR_doublets.txt"), col.names = TRUE, row.names = TRUE, sep = "\t")

    }
}

#################################################################################
