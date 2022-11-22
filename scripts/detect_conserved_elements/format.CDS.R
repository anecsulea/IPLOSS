## original script by Alexandre Laverré & Anamaria Necsulea

######################################################################

pathAnnot="../../data/ensembl_annotations/" ## define path
prefixAnnot="AllTranscripts_Ensembl103" ## file name for the annotations

spnames=c("Gallus_gallus", "Anas_platyrhynchos_platyrhynchos")
names(spnames)=c("Chicken", "Duck")

######################################################################

for(sp in c("Chicken", "Duck")){
    name=spnames[sp]

    ens=read.table(paste(pathAnnot,sp, "/", prefixAnnot,".gtf",sep=""), h=F, sep="\t", stringsAsFactors=F, quote="\"")
    ens=ens[grep("transcript_biotype protein_coding", ens$V9),]
    cds=ens[which(ens$V3=="CDS"),]

    cds$V1=paste(name, cds$V1, sep=".") ## change chromosome names so that they match those from the MAF file

    for (chr in c(as.character(1:29),"Z","W")){
        print(chr)
        nb=length(which(cds$V1 == paste(name, chr, sep=".")))

        if(nb>0){
            this.cds=cds[which(cds$V1 == paste(name, chr, sep=".")),]
            write.table(this.cds, file=paste(pathAnnot, sp, "/by_chr/", prefixAnnot,".CDS_",chr,".gtf",sep=""), row.names=F, col.names=F, sep="\t", quote=F)
        }
    }
}


######################################################################
