#######################################################################################

## adapted from original script by Menghan Wang
## follow the protocol of this two pages:
## https://medium.com/@lokrajvet/working-with-the-genome-of-non-model-organism-in-r-bioconductor-8c6eec253c8a
## https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf

#######################################################################################

library(ArchR)
library(BSgenome)
library(GenomeInfoDb)
library(ensembldb)

#######################################################################################

set.seed(1234)
addArchRThreads(threads = 1)

#######################################################################################

abbr=c("Ggallus", "Aplatyrhynchos")
names(abbr)=c("Chicken", "Duck")

pathGenomeIndexes="../../results/snATACSeq_indexes/"

#######################################################################################

for(sp in c("Chicken", "Duck")){

  this.abbr = abbr[sp]
  
  my_file = read.dcf(paste(pathGenomeIndexes, sp, "/seed_file.txt", sep=""), fields = NULL, all = FALSE, keep.white = TRUE)
  write.dcf(my_file, file = paste(pathGenomeIndexes, sp, "/seed.dcf", sep=""), append = FALSE, useBytes = FALSE, indent = 0.1 * getOption("width"), width = 0.9 * getOption("width"), keep.white = NULL)
            
  forgeBSgenomeDataPkg(paste(pathGenomeIndexes, sp, "/seed.dcf",sep=""))
  
  system(paste("R CMD build BSgenome.", this.abbr, sep=""))
  system(paste("R CMD check BSgenome.", this.abbr,"_1.0.0.tar.gz", sep=""))

  system(paste("R CMD INSTALL BSgenome.", this.abbr,"_1.0.0.tar.gz", sep=""))
  
}

#######################################################################################
