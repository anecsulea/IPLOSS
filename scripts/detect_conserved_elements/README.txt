######################################################################

## Step 0: prepare the input alignments for your set of species of interest.

You can use mafSpeciesSubset from UCSC utilities
(http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) to extract a
subset of species from a maf alignment. 

Usage:  mafSpeciesSubset in.maf species.lst out.maf

######################################################################

Step 1: extract the phylogenetic tree for the same set of species.

- select.IP.species.R 

You can start
with the newick tree for the full set of species
(e.g. https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/multiz77way/galGal6.77way.nh),
and then use the keep.tip command in the ape library in R to keep only
the species that you're interested in (or drop.tip to remove the other
ones if it's easier).

######################################################################

## Step 2: format the CDS annotations to match the alignemnts. 

- format.CDS.R

 Adapt this script to extract the CDS from a reference GTF file, for example from
Ensembl, and make sure that the species names in the alignments appear
with the same convention in the
GTF files. In the MAF files from UCSC, the species names and
chromosomes are named as in this example:

a score=-5.000000
s galGal6.chrM     0 7 + 16775 aatttta
s melGal5.chrM 15553 7 + 16719 aataatt
i melGal5.chrM N 0 C 0

In Ensembl, the chromosomes are named 1,2,3, Z, W, MT... You need to
rename them as follows:

   1 -> galGal6.chr1 
   2 -> galGal6.chr2
   ... 
   MT -> galGal6.chrM (maybe just discard this one).

It will be useful for later to output the CDS coordinates in separate files, one
for each chromosome.

######################################################################

## Step 3: prepare the alignments for 4-fold degenerate sites (neutrally evolving regions)

- extract.4fold.degenerate.sites.sh

This script uses msa_view (which comes with the phastCons package) to
extract the alignments corresponding to the 4-fold degenerate
sites. We do this separately for each chromosome. There are two
msa_view commands, one to extract the codons which contain 4-fold
degenerate sites, and then the 4-fold sites themselves (third codon
positions). 

######################################################################

## Step 4: estimate the rates of evolution for 4-fold degenerate sites with phyloFit

- run.phyloFit.sh

phyloFit (from the phastCons package) needs the phylogenetic tree of the species of interest
(generated at step 1) and the alignments for the 4-fold degenerate
sites (generated at step 3). We do this separately for each chromosome.

######################################################################

## Step 5: combine the phyloFit evolutionary models for chromosome
   types (autosomal macro-chromosomes, autosomal micro-chromosomes, sex chromosomes)

- combine.phyloFit.sh

We use phyloBoot from the phastCons package to do this.
You can check out http://www.ensembl.org/Gallus_gallus/Location/Genome
for an overview of the karyotype. The boundary between macro and
micro-chromosomes is not obvious. Based on some definitions from the
literature, we decided to call macro-chromosomes those autosomes that
are larger than 25Mb (so chromosomes 1 to 8); microchromosomes are
chromosomes 9 to 39 for chicken. You can combine the Z and the W
chromosomes for the sex chromosome model.

######################################################################

## Step 6: run phastCons to identify conserved elements and obtain
   per-base conservation scores

- run.phastCons.sh 

phastCons will need as input the evolutionary models (obtained at step
5) for neutrally evolving sites and the MAF files. There are three
parameters that you can play with: --target-coverage ,
--expected-length and --rho.
Here we are estimating the rho values, so we just give a starting
value for the estimation. We ask for two outputs: the per-base
conservation score in wig format (by using the --score parameter), and
the coordinates of the conserved elements (--most-conserved).

NB: there are other ways to run phastCons. Here we estimate the
evolutionary model for neutrally evolving positions and we ask
phastCons to estimate the model for the conserved positions, through
the rho parameter (ratio of the rate of evolution between conserved
and non-conserved elements). Another option would be to also provide
phastCons with the evolutionary model for conserved positions (e.g.,
non-synonymous positions in CDS), that would mean that the conserved
elements we obtain at the end would be as conserved as coding
exons. See also http://compgen.cshl.edu/phast/phastCons-HOWTO.html . 

Check out section 3.3 in the phastCons doc for more advice on how to tune
the parameters. 

######################################################################

## Step 7: compute average phastCons scores for the most conserved
   elements (optional)

- compute.average.scores.sh

This uses the home-made compute.average.scores.pl script, to compute
average scores for the elements obtained at step 6 (or for any other
elements in bed format). This script uses the .wig file containing the
per-base conservation score obtained at step 6. I don't know what the
4th column in the bed files obtained at step 6 corresponds to, but
it's not an average score, maybe the sum of the scores or the covered fraction?

######################################################################
