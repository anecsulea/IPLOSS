#####################################################################

library(Seurat)

#####################################################################

pathResults="../../results/snRNASeq_analysis/"
annot="EnsemblStringTie"

#####################################################################

samples <- list()
samples[["Chicken"]] <- "Gg_W_GT_male_33_scRNA"
samples[["Duck"]] <- c("Ap_W_GT_male_33_scRNA", "Ap_W_GT_female_33_scRNA")

#####################################################################

for(sp in names(samples)){
    for(sample in samples[[sp]]){

        load(paste(pathResults, sp, "/",annot, "/",sample, "/seurat.object.after.clustering.RData",sep=""))

        so <- filtered.seurat

        ## prepare data for CellPhoneDB
        count_raw <- GetAssayData(object=so, slot="counts")
        count_raw <- as.matrix(count_raw)

        count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
        write.table(count_norm, paste(pathResults, "data_for_cellphonedb/cellphonedb_count_",sp,"_",sample,"_",annot,".txt",sep=""), sep="\t", quote=F)

        ## generating meta file
        all_meta <- so[[]]
        meta_data <- cbind(rownames(all_meta), all_meta[,c("seurat_clusters")])
        write.table(meta_data,  paste(pathResults, "data_for_cellphonedb/cellphonedb_meta_",sp,"_",sample,"_", annot, ".txt",sep=""), sep="\t", quote=F, row.names=F)
    }
}

#####################################################################
