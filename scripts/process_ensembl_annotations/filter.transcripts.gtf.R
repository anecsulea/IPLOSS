##############################################################################

release=103

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

maxlen=2500000

##############################################################################

for(sp in c("Chicken", "Duck")){
  pathAnnot=paste("../../data/ensembl_annotations/",sp, sep="")

  ## readthrough transcripts

  rt=read.table(paste(pathAnnot,  "/ReadthroughTranscripts_Ensembl",release, ".txt", sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

#################################

  geneinfo=read.table(paste(pathAnnot,  "/GeneInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")

#################################

  txinfo=read.table(paste(pathAnnot,  "/TranscriptInfo_Ensembl",release,".txt",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="\"")
  colnames(txinfo)[1]="gene"
  colnames(txinfo)[2]="tx"

#################################

  monoexonic=txinfo$tx[which(txinfo$biotype%in%c("snoRNA", "miRNA", "snRNA"))]
  print(paste(length(monoexonic), "monoexonic transcripts (miRNA, snRNA, snoRNA)"))

#################################

  txinfo$Length=txinfo[,"seq_region_end"]-txinfo[,"seq_region_start"]+1
  oklen=txinfo$tx[which(txinfo$Length<=maxlen)]

#################################

  ## for human, some genes have transcripts on multiple strands, we remove them - trans-splicing events, not sure StringTie can handle them

  nbstrands=tapply(txinfo$seq_region_strand, as.factor(txinfo$gene), function(x) length(unique(x)))
  weirdgenes=names(nbstrands)[which(nbstrands!=1)]

  if(length(weirdgenes)>0){
    print(paste(length(weirdgenes), "genes with more than one strand, likely trans-splicing"))
    print(weirdgenes)
  }

#################################

  print(paste(length(which(txinfo$tx%in%rt$TranscriptID)), " read-through transcripts"))

  ok=which((txinfo$tx%in%oklen) & (!(txinfo$tx%in%rt$TranscriptID)) & (!txinfo$gene%in%weirdgenes) & (!txinfo$tx%in%monoexonic))

  selected=txinfo[ok,2]

  print("We remove read-through transcripts, trans-splicing events, genes that are too long (>2.5Mb) and monoexonic transcripts - to be added later.")
  print(paste(length(selected), "selected transcripts"))

#################################

  print("reading gtf")

  gtf=read.table(paste(pathAnnot,  "/AllTranscripts_Ensembl",release,".gtf",sep=""), h=T, stringsAsFactors=F, sep="\t", quote="")

  colnames(gtf)=c("chr", "source", "type", "start", "end", "score1", "strand", "score2", "info")

  gtf=gtf[which(gtf$type=="exon"),]

  print("done")

##################################

  gtf.info=lapply(gtf$info, function(x) unlist(strsplit(x, split=";")))
  txid=unlist(lapply(gtf.info, function(x) grep("transcript_id", x, value=T)))
  txid=unlist(lapply(txid, function(x) unlist(strsplit(x,split="\""))[2]))

  print(paste(length(unique(txid)), "unique transcripts in GTF file"))

###################################

  gtf.selected=gtf[which(txid%in%selected),]

  write.table(gtf.selected, file=paste(pathAnnot, "/FilteredTranscripts_Ensembl",release,".gtf",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

###################################

  gtf.monoex=gtf[which(txid%in%monoexonic),]

  write.table(gtf.monoex, file=paste(pathAnnot, "/monoexonicTranscripts_Ensembl",release,".gtf",sep=""), row.names=F, col.names=F, sep="\t", quote=F)

}

##############################################################################
