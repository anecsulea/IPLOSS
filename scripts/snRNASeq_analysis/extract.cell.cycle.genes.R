## cell cycle markers from Scialdone et al., 2015
## orthologues for chicken extracted by Menghan Wang, University of Basel

######################################################################

pathCellCycle="../../data/single_cell_data/cell_cycle_genes/"
pathHomology="../../data/ensembl_homology/"

release="103"

######################################################################

load(paste(pathCellCycle, "chicken_cycle_markers.Rda",sep=""))

######################################################################

enshom=read.table(paste(pathHomology, "HomologousGenes_Chicken_Duck_Ensembl",release,".txt",sep=""),h=T,stringsAsFactors=F, sep="\t", quote="\"")

enshom=enshom[which(enshom$Duck.homology.type=="ortholog_one2one"),]
rownames(enshom)=enshom$Gene.stable.ID

######################################################################

ap.pairs=list()

for(phase in names(gg.pairs)){
    this.pairs=gg.pairs[[phase]]

    ortho.first=enshom[this.pairs$first,"Duck.gene.stable.ID"]
    ortho.second=enshom[this.pairs$second,"Duck.gene.stable.ID"]

    ok=which((!is.na(ortho.first)) & !is.na(ortho.second))
    this.ap.pairs=data.frame("first"=ortho.first[ok], "second"=ortho.second[ok], stringsAsFactors=F)

    ap.pairs[[phase]]=this.ap.pairs

    print(paste(phase, nrow(this.pairs), "pairs originally", nrow(this.ap.pairs), "ortho pairs"))
}

######################################################################

save(ap.pairs, file=paste(pathCellCycle, "duck_cycle_markers.RDa",sep=""))

######################################################################
