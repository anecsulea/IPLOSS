########################################################################

set.seed(19)

options(stringsAsFactors=F)

library(tximport)

########################################################################

pathExpression="../../results/gene_expression_estimation/"

ensrelease=103

types=c(paste("AllTranscripts_Ensembl", ensrelease, sep=""), "EnsemblStringTie")
splist=c("Chicken", "Duck")

########################################################################

source("normalization.R")

########################################################################

for(sp in splist){

  for(type in types){
    print(paste(sp, type))

    if(file.exists(paste(pathExpression,  sp, "/", type,"/", "AllSamples_KallistoEffectiveLength.txt",sep=""))){
      print("already done")
    } else{
      ## sample list

      samples=system(paste("ls ", pathExpression, sp, "/", type," | grep -v txt | grep GT", sep=""), intern=T)

      if(length(samples)>0){

          sinfo=read.table(paste(pathExpression, sp, "/", type, "/",samples[1],"/abundance.tsv",sep=""), h=T, stringsAsFactors=F)
          geneid=unlist(lapply(sinfo$target_id, function(x) unlist(strsplit(x, split=":"))[1]))

          tx2gene=data.frame("txid"=sinfo$target_id, "geneid"=geneid)

          ##############################################################################

          files=paste(pathExpression, sp, "/", type,"/", samples, "/abundance.tsv", sep="")
          names(files)=samples

          txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)

          ########################################################################

          read.counts=as.data.frame(txi.kallisto$counts)

          tpm=as.data.frame(txi.kallisto$abundance)

          efflen=as.data.frame(txi.kallisto$length)

          ########################################################################

          ## normalization

          norm.data=normalization(tpm)
          tpm.norm=norm.data[["expdata.norm"]]
          rownames(tpm.norm)=rownames(tpm)

          hk.genes=norm.data[["hk.genes"]]

          ########################################################################

          ## add gene id as a column

          tpm$GeneID=rownames(tpm)
          tpm.norm$GeneID=rownames(tpm.norm)
          read.counts$GeneID=rownames(read.counts)
          efflen$GeneID=rownames(efflen)

          tpm=tpm[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
          tpm.norm=tpm.norm[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
          read.counts=read.counts[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]
          efflen=efflen[,c("GeneID", setdiff(colnames(tpm), "GeneID"))]

        ########################################################################

        ## write output

          write.table(efflen, file=paste(pathExpression,  sp, "/", type,"/", "AllSamples_KallistoEffectiveLength.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

          write.table(read.counts, file=paste(pathExpression, sp, "/", type,"/",  "AllSamples_KallistoEstimatedCounts.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

          write.table(tpm, file=paste(pathExpression,   sp, "/", type,"/", "AllSamples_KallistoRawTPM.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

          write.table(tpm.norm, file=paste(pathExpression,  sp, "/", type,"/",  "AllSamples_KallistoNormalizedTPM.txt",sep=""), row.names=F, col.names=T, quote=F, sep="\t")

       ########################################################################
      }
    }
  }
}

########################################################################
