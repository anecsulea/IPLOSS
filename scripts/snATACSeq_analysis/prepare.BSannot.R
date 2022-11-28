
##########################
#Creating a Custom ArchRGenome
##########################

#https://www.archrproject.com/bookdown/getting-set-up.html
library(BSgenome.Ggallus.ensembl.galGal6) #library(BSgenome.Ggallus.UCSC.galGal6)
# in order to keep chrZ and chrW, first change filter=TRUE to filter=FALSE
genomeAnnotation = createGenomeAnnotation(genome = "BSgenome.Ggallus.ensembl.galGal6",filter=FALSE)
# then manually retain autosomes and sex chromosomes
genomeAnnotation$chromSizes=genomeAnnotation$chromSizes[2:35]
unique(seqnames(genomeAnnotation$chromSizes));seqlengths(BSgenome.Ggallus.ensembl.galGal6)

##gtf need to have gene_name
#awk -F $'\t' '$9 ~ /gene_name/ {print $0}' Gg6_extended_060420.gtf > Gg6_extended_060420.filter.gtf

anno_path <- "/scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/"
db <- ensDbFromGtf(gtf=paste0(anno_path,"Gg6_extended_060420.filter.gtf"),
                   organism = "Gallus_gallus", genomeVersion = "GRCg6a", version = 97)
edb <- EnsDb(db); rm(db)
edb = addFilter(edb, filter = GeneBiotypeFilter('protein_coding',"=="))
# change the settings here as well, to keep consistent
edb = addFilter(edb, filter = SeqNameFilter(c('AADN05001525.1','KZ626819.1','KZ626826.1','KZ626830.1','MT',
                                              'KZ626833.1','KZ626834.1','KZ626835.1','KZ626836.1','KZ626839.1'),"!=")) #'Z','W',
gene.ranges <- genes(edb); length(gene.ranges$gene_id)
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges),
                      ranges = IRanges(start = start(gene.ranges), width = 2),
                      strand = strand(gene.ranges)
)
exon.ranges = exons(edb)
geneAnnotation = createGeneAnnotation(
  TSS = tss.ranges,
  exons = exon.ranges,
  genes = gene.ranges
)
unique(seqnames(geneAnnotation$genes))

save_path = "/scicore/home/tschoppp/wang0007/genitalia/"
saveRDS(genomeAnnotation, paste0(save_path,"genomeAnnotation_withZW.rds"))
saveRDS(geneAnnotation, paste0(save_path,"geneAnnotation_withZW.rds"))



