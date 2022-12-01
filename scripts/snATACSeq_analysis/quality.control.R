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

pathResults="../../results/snATACSeq_analysis/"
annot="Ensembl"

samples=list("Chicken"=c("Gg_W_GT_male_33_scATAC"), "Duck"=c("Ap_W_GT_male_33_scATAC"))

#################################################################################

for(sp in c("Chicken", "Duck")){

    load(paste(pathResults, sp, "/",annot,"/geneAnnotation.RData",sep=""))

    for(sample in samples[[sp]]){

        path.res = paste(pathResults, sp, "/", annot, "/",sample, "/outs/", sep="")

        ## input to Signac

        counts = Read10X_h5(filename = paste0(path.res, "filtered_peak_bc_matrix.h5"));
        metadata = read.csv(file = paste0(path.res, "singlecell.csv"),header = TRUE,row.names = 1)
        chrom_assay = CreateChromatinAssay(counts = counts, sep = c(":", "-"), min.cells = 10, min.features = 200)

        ## create Seurat object

        mtx <- CreateSeuratObject(
            counts = chrom_assay,
            assay = "peaks",
            meta.data = metadata
        )

        ## fragment object

        frags <- CreateFragmentObject(path = paste0(path.res,'fragments.tsv.gz'), cells = colnames(x = mtx), validate.fragments = FALSE)
        Fragments(mtx[["peaks"]]) <- frags

        ## divide genome into bins, add them to count matrix

        if(sp=="Chicken"){
            library(BSgenome.Ggallus)
            genome = seqlengths(BSgenome.Ggallus)
        }

         if(sp=="Duck"){
             library(BSgenome.Aplatyrhynchos)
             genome = seqlengths(BSgenome.Aplatyrhynchos)
         }

        genome=genome[names(genome) %in% c(1:33,"Z")]

        binmatrix = GenomeBinMatrix(
            fragments = Fragments(mtx),
            genome = genome,
            cells = colnames(mtx),
            binsize = 5000, verbose = F
        )

        mtx[['bins']] <- CreateAssayObject(counts = binmatrix)

        ## get previously saved gene annotation

        gene.ranges = geneAnnotation[["genes"]]
        annotations = gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]

        ## add gene annotation to matrix
        Annotation(mtx) <- annotations

        ## nucleosome signal
        mtx <- NucleosomeSignal(object = mtx,n = ncol(mtx) * 10000)

        ## TSS enrichment

        mtx <- TSSEnrichment(object = mtx, fast = FALSE)
        mtx$high.tss <- ifelse(mtx$TSS.enrichment > 2, 'High', 'Low')
        TSSPlot(mtx, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

        ## percentage of reads in peaks
        mtx$pct_reads_in_peaks <- mtx$peak_region_fragments / mtx$passed_filters * 100

        ## filter matrix
        ## peak_region_fragments between 500 and 100000
        ## more than 5% reads in peaks
        ## nucleosome signal < 4
        ## TSS enrichment > 1
        mtx.filtered <- subset(mtx, subset = peak_region_fragments > 500 &
                   peak_region_fragments < 100000 &
                   pct_reads_in_peaks > 5 &
                   nucleosome_signal < 4 & TSS.enrichment > 1)

        ## filter also on doublets identified by ArchR

        doublets=read.table(paste0(pathResults, sp, "/", annot, "/",sample,"/ArchR_doublets.txt"), header = TRUE, row.names = 1, sep = "\t", stringsAsFactors=F)

        doublets=doublets[,1]
        doublets=unlist(lapply(doublets, function(x) unlist(strsplit(x, split="\\#"))[2]))

        mtx.filtered@meta.data$doublet = "F"
        mtx.filtered@meta.data[rownames(mtx.filtered@meta.data) %in% doublets,"doublet"]="T"

        save(list=c("mtx", "mtx.filtered"), file=paste(pathResults, sp, "/", annot, "/",sample,"/matrix.afterQC.RData",sep=""))
    }

}

#################################################################################
