######################################################################

pathStringTie="../../results/stringtie_assembly/"
pathEnsembl="../../data/ensembl_annotations/"

release=103

options(scipen=999) ## remove scientific notation ## options(scipen=0) to get it back

######################################################################

for(sp in c("Chicken", "Duck")){

  print(sp)

  st=read.table(paste(pathStringTie, sp, "/ExonBlocks_combined_annotations_StringTie_Ensembl.txt", sep=""), h=F, stringsAsFactors=F, sep="\t", quote="")

  st=st[,c(3,4,5,6)]
  colnames(st)=c("chr", "start", "end", "strand")

  st$id=paste(st$chr, st$start, st$end, st$strand, sep=",")

  print(dim(st)[1])

  st$score=rep("1000", dim(st)[1])
  st$start=st$start-1 ## bed format

  dupli=which(duplicated(st$id))

  if(length(dupli)>0){
    st=st[-dupli,]
    print(paste("removed", length(dupli), "duplicated lines"))
  }

  st$outstrand=rep(NA, dim(st)[1])
  st$outstrand[which(st$strand%in%c("1", "+"))]="+"
  st$outstrand[which(st$strand%in%c("-1", "-"))]="-"

  print(paste(length(which(is.na(st$outstrand))), "NA values for strand"))

  write.table(st[,c("chr", "start", "end", "id", "score", "outstrand")], file=paste(pathStringTie, sp, "/ExonBlocks_combined_annotations_StringTie_Ensembl.bed", sep=""), row.names=F, col.names=F, sep="\t", quote=F)
}

######################################################################
