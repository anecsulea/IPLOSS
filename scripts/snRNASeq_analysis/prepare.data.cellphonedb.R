#####################################################################

library(Seurat)

#####################################################################

pathResults="../../results/snRNASeq_analysis/"
pathHomology="../../data/ensembl_homology/"

#####################################################################

samples <- list()
samples[["Chicken"]] <- "Gg_W_GT_male_33_scRNA"
samples[["Duck"]] <- c("Ap_W_GT_male_33_scRNA", "Ap_W_GT_female_33_scRNA")

#####################################################################

for(sp in names(samples)){
    enshom  <-  read.table(paste(pathHomology, "HomologousGenes_Human_",sp,"_Ensembl103.txt",sep=""),h=T, stringsAsFactors=F, sep="\t", quote="\"")
    enshom  <-  enshom[which(enshom[,3]=="ortholog_one2one"),]
    rownames(enshom)  <-  enshom[,2]

    for(sample in samples[[sp]]){

        for(annot in c("Ensembl", "EnsemblStringTie")){

            system(paste("mkdir -p ", pathResults, sp, "/",annot,"/",sample, "/cellphonedb",sep=""))

            load(paste(pathResults, sp, "/",annot, "/",sample, "/seurat.object.after.clustering.RData",sep=""))

            so <- filtered.seurat

            ## prepare data for CellPhoneDB
            count_raw <- GetAssayData(object=so, slot="counts")
            count_raw <- as.matrix(count_raw)

            ## take only ortho genes
            count_raw <-  count_raw[which(rownames(count_raw)%in%rownames(enshom)),]
            rownames(count_raw) <-  enshom[rownames(count_raw),1]

            count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
            count_norm  <-  as.data.frame(count_norm)
            count_norm$Gene  <-  rownames(count_raw)
            count_norm  <- count_norm[,c("Gene", colnames(count_raw))]

            write.table(count_norm, paste(pathResults, sp, "/",annot,"/",sample, "/cellphonedb/cellphonedb_count.txt",sep=""), sep="\t", quote=F, row.names=F)

            ## generating meta file
            all_meta <- so[[]]
            meta_data <- cbind(rownames(all_meta), paste("ctype",all_meta[,c("seurat_clusters")],sep=""))
            colnames(meta_data) <- c("Cell", "cell_type")

            write.table(meta_data,  paste(pathResults,sp, "/", annot, "/",sample, "/cellphonedb/cellphonedb_meta.txt",sep=""), sep="\t", quote=F, row.names=F)
        }
    }
}

#####################################################################
