######################################################################

pathStringTie="../../results/stringtie_assembly/"
pathExpression="../../results/expression_estimation/"

mode="reference_none_selected_samples"

######################################################################

for(sp in c("Chicken", "Duck")){
  tpm.ens=read.table(paste(pathExpression, sp, "/AllTranscripts_Ensembl103/GTSamples_KallistoRawTPM.txt",sep=""),h=T, stringsAsFactors=F)
  rownames(tpm.ens)=tpm.ens$GeneID
  tpm.ens=tpm.ens[,which(colnames(tpm.ens)!="GeneID")]

  stage=unlist(lapply(colnames(tpm.ens), function(x) unlist(strsplit(x, split="_"))[4]))
  tpm.ens=tpm.ens[,which(stage!="31" & stage!="32")]

  sumtpm.ens=apply(tpm.ens, 1, sum)
  maxtpm.ens=apply(tpm.ens, 1, max)

  se=read.table(paste(pathStringTie, sp, "/", mode,"/assembled_transcripts_vs_Ensembl.txt",sep=""), h=T, stringsAsFactors=F)
  es=read.table(paste(pathStringTie, sp, "/", mode,"/Ensembl_vs_assembled_transcripts.txt",sep=""), h=T, stringsAsFactors=F)

  se$FractionOverlap1=se$LengthOverlap/se$LengthTranscript1
  se$FractionOverlap2=se$LengthOverlap/se$LengthTranscript2

  es$FractionOverlap1=es$LengthOverlap/es$LengthTranscript1
  es$FractionOverlap2=es$LengthOverlap/es$LengthTranscript2

  se$PartialMatch=(se$FractionOverlap1>=0.5 & se$FractionOverlap2>=0.5)
  se$GoodMatch=(se$FractionOverlap1>=0.9 & se$FractionOverlap2>=0.9)

  es$PartialMatch=(es$FractionOverlap1>=0.5 & es$FractionOverlap2>=0.5)
  es$GoodMatch=(es$FractionOverlap1>=0.9 & es$FractionOverlap2>=0.9)

  maxov.se.tx=tapply(se$FractionOverlap1, as.factor(se$TranscriptID1), max)
  maxov.es.tx=tapply(es$FractionOverlap1, as.factor(es$TranscriptID1), max)

  maxov.se.gene=tapply(se$FractionOverlap1, as.factor(se$GeneID1), max)
  maxov.es.gene=tapply(es$FractionOverlap1, as.factor(es$GeneID1), max)

  sumtpm.ens=sumtpm.ens[names(maxov.es.gene)]

  tpm.classes=cut(sumtpm.ens, breaks=c(seq(from=0, to=50, by=5), max(sumtpm.ens)), include.lowest=T)

  boxplot(maxov.es.gene~tpm.classes)

  weird=names(maxov.es.gene)[which(maxov.es.gene==0 & sumtpm.ens>100)]

}
######################################################################
