## original script by Menghan Wang, University of Basel

###########################################################################

pathResults="../../results/snRNASeq_analysis/"
pathAnnot="../../data/ensembl_annotations/"

splist=c("Chicken", "Duck")
annot=c("Ensembl", "EnsemblStringTie")
samples=list("Chicken"=c("Gg_W_GT_male_33_scRNA"), "Duck"=c("Ap_W_GT_female_33_scRNA",  "Ap_W_GT_male_33_scRNA"))
sample.names=c("gg_male", "ap_female", "ap_male")
names(sample.names)=c("Gg_W_GT_male_33_scRNA", "Ap_W_GT_female_33_scRNA",  "Ap_W_GT_male_33_scRNA")

###########################################################################

## load libraries

library(Seurat)

###########################################################################

for(sp in splist){
    ## read gene info and select mitochondrial genes
    gene.info=read.table(paste(pathAnnot, sp, "/GeneInfo_Ensembl103.txt",sep=""), header = TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
    mt.genes=gene.info$stable_id[which(gene.info$name=="MT")]

    for(sample in samples[[sp]]){
        for(annot in c("Ensembl", "EnsemblStringTie")){

            print(paste(sp, sample, annot))

            short.name=sample.names[sample]

            ## save all plots here
            pdf(file=paste(pathResults, sp, "/",annot, "/",sample, "/QualityControls.pdf",sep=""))

            ## create Seurat object

            all.data=Read10X(data.dir = paste0(pathResults,sp, "/",annot, "/", sample, "/outs/filtered_feature_bc_matrix/"), gene.column = 1)

            ## initial statistics
            initial.stats = data.frame(cells=colnames(all.data),UMIs=Matrix::colSums(all.data), genes = apply(all.data,2, function(x) length(which(x>0))))

            par(mfrow=c(2,2))
            par(oma=c(0,0,2,0))
            hist(initial.stats$UMIs, xlab="UMIs/cell", ylab="Frequency", breaks=30, main="")
            hist(initial.stats$genes/initial.stats$UMIs, xlab="UMIs/gene", ylab="Frequency", breaks=30, main="")
            hist(initial.stats$genes, xlab="Number of genes per cell", ylab="Frequency", breaks=30, main="")
            plot(initial.stats$UMIs, initial.stats$genes, pch=20, main="", xlab="Number of UMIs", ylab="Number of genes")

            mtext("before filtering", outer=T, line=1, side=3, font=2)
            par(mfrow=c(1,1))

            ctable=all.data ## we will filter out cells from this object, all.data contains all cells
            print(paste(ncol(ctable), "cells originally"))

            ## plot the distribution of UMIs per cell
            ## stat.umi = data.frame(cells=colnames(all.data),UMIs=Matrix::colSums(all.data))
            ## hist(stat.umi$UMIs, xlab="UMIs/Cell", ylab="Frequency",main=paste0(short.name," Raw UMIs/Cell Frequency"), breaks=30)
            ## abline(v=mean(stat.umi$UMIs))

            ## We first filter out the cells that have more than 4 times the mean of the UMI counts, and less than 20% the mean.
            ## These represent the probable doublets and poor cells.

            ## Remove those cells that have more than 4 times the mean of UMI count
            ctable = ctable[,Matrix::colSums(ctable)<(4*mean(Matrix::colSums(ctable)))]
            print(paste(ncol(ctable), "cells after removing probable doublets"))

            ## # We remove the cells that have less than 20% of the UMIs median because they're likely false positives
            ctable = ctable[,Matrix::colSums(ctable)>(0.2*median(Matrix::colSums(ctable)))]
            print(paste(ncol(ctable), "cells after removing poor quality cells (false positives)"))

            ## The second filter is for mitochondrial UMIs, we get rid of cells that show that their fraction of UMIs from MT origin is higher than:
            ## The median + three times the MAD (median absolute deviation). In this case `r round(median(percent.mt) + (mad(percent.mt) * 3),3)`

            is.mt = (rownames(ctable) %in% mt.genes)
            prop.mt = Matrix::colSums(ctable[is.mt, ])/Matrix::colSums(ctable)
            nb.UMIs=Matrix::colSums(ctable)

            plot(prop.mt, nb.UMIs, xlab="MT proportion", ylab="UMI counts", pch=20)

            my.threshold = median(prop.mt) + (mad(prop.mt) * 3)
            filter.mtprop = rep(1,length(prop.mt)) ## 1 if we keep cell, 0 if not

            ## Original filter: we remove genes that have a MT proportion higher than the threshold, among those genes that are in the bottom 50% in terms of number of UMIs
            ## filter.mtprop[which(prop.mt> my.threshold & Matrix::colSums(ctable) < median(Matrix::colSums(ctable)) )] = 0

            ## Slightly modified filter: we remove genes that have a MT proportion higher than the threshold, irrespective of the number of UMIs that they have
            filter.mtprop[which(prop.mt> my.threshold)] = 0

            points(prop.mt[filter.mtprop==0], nb.UMIs[filter.mtprop==0], col="red", pch=20)

            ctable=ctable[,which(filter.mtprop==1)]
            print(paste(ncol(ctable), "cells after removing cells with a high proportion of MT reads (>", my.threshold,")"))

            ## The next filter is for cells that have an unusual number of genes detected per UMI

            ## Number of genes detected as expressed in each cell
            gct = data.frame(cell = colnames(ctable), genes = apply(ctable,2, function(x) length(which(x>0))))

            ## Column for UMI counts
            gct = cbind(gct, UMIs=(Matrix::colSums(ctable)))

            ## plot(gct$UMIs, gct$genes, xlab="nb UMIs", ylab="nb detected genes", pch=20)

            filter.geneumis=rep(1, nrow(gct))

            ## We remove cells for which the ratio of genes to UMIs is below 0.15
            filter.geneumis[which((gct$genes/gct$UMIs)<0.15)]=0

            ctable=ctable[,which(filter.geneumis==1)]
            print(paste(ncol(ctable), "cells after removing cells for which the proportion of genes to UMIs is below 0.15"))

            ## final statistics plot

            final.stats = data.frame(cells=colnames(ctable),UMIs=Matrix::colSums(ctable), genes = apply(ctable,2, function(x) length(which(x>0))))

            par(mfrow=c(2,2))
            par(oma=c(0,0,2,0))
            hist(final.stats$UMIs, xlab="UMIs/cell", ylab="Frequency", breaks=30, main="")
            hist(final.stats$genes/final.stats$UMIs, xlab="UMIs/gene", ylab="Frequency", breaks=30, main="")
            hist(final.stats$genes, xlab="Number of genes per cell", ylab="Frequency", breaks=30, main="")
            plot(final.stats$UMIs, final.stats$genes, pch=20, main="", xlab="Number of UMIs", ylab="Number of genes")
            mtext("after filtering", outer=T, line=1, side=3, font=2)

            filtered.data = all.data[,colnames(ctable)]

            ## save objects

            filtered.seurat=CreateSeuratObject(filtered.data,project = short.name)
            filtered.seurat@meta.data$percent.mt =  PercentageFeatureSet(filtered.seurat, features = intersect(mt.genes, rownames(all.data)))

            save(filtered.seurat, file=paste(pathResults, sp, "/",annot, "/",sample, "/seurat.object.afterQC.RData",sep=""))

            ## close pdf
            dev.off()
        }
    }
}

###########################################################################
