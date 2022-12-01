#################################################################################

## adapted from original script by Menghan Wang
## https://www.archrproject.com/bookdown/getting-set-up.html

#################################################################################

library(ArchR)
library(BSgenome)
library(GenomeInfoDb)
library(ensembldb)

#################################################################################

pathEnsembl="../../data/ensembl_annotations/"
pathEnsemblStringTie="../../results/stringtie_assembly/"
pathResults="../../results/snATACSeq_analysis/"

#################################################################################

for(sp in c("Chicken", "Duck")){

    if(sp=="Chicken"){
        library(BSgenome.Ggallus)

        genomeAnnotation = createGenomeAnnotation(genome = "BSgenome.Ggallus",filter=FALSE)

        ## we remove mitochondria, W chromosome, unassembled contigs
        genomeAnnotation$chromSizes=genomeAnnotation$chromSizes[c(1:32, 34)]

        organism="Gallus_gallus"
    }

    if(sp=="Duck"){
        library(BSgenome.Aplatyrhynchos)

        genomeAnnotation = createGenomeAnnotation(genome = "BSgenome.Aplatyrhynchos",filter=FALSE)

        organism="Anas_platyrhynchos_platyrhynchos"
    }

    chrlist = as.character(unique(seqnames(genomeAnnotation$chromSizes)))

    path.annot = paste(pathEnsembl, sp, "/AllTranscripts_Ensembl103.gtf",sep="")

    if(sp=="Chicken"){
        edb = ensDbFromGtf(gtf=path.annot, organism = organism, genomeVersion = "GRCg6a", version = 103)
    }

    if(sp=="Duck"){
        edb = ensDbFromGtf(gtf=path.annot, organism = organism, genomeVersion = "CAU_duck1.0", version = 103)
    }

    edb = EnsDb(edb)

    ## keep protein-coding genes
    edb = addFilter(edb, filter = GeneBiotypeFilter('protein_coding',"=="))

    ## same chromosomes as above
    edb = addFilter(edb, filter = SeqNameFilter(chrlist,"=="))

    gene.ranges = genes(edb)

    gene.strand=strand(gene.ranges)
    gene.start=start(gene.ranges)
    gene.end=end(gene.ranges)

    tss.pos=rep(NA, length(gene.start))
    tss.pos[which(gene.strand=="+")]=gene.start[which(gene.strand=="+")]
    tss.pos[which(gene.strand=="-")]=gene.end[which(gene.strand=="-")]

    tss.ranges = GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = tss.pos, width = 2), strand = strand(gene.ranges))

    exon.ranges = exons(edb)

    geneAnnotation = createGeneAnnotation(TSS = tss.ranges, exons = exon.ranges, genes = gene.ranges)

    print(unique(seqnames(geneAnnotation$genes)))

    save(genomeAnnotation, file=paste(pathResults, sp, "/Ensembl/genomeAnnotation.RData",sep=""))
    save(geneAnnotation, file=paste(pathResults, sp, "/Ensembl/geneAnnotation.RData",sep=""))
}

#################################################################################
