#!/bin/bash

export sp=$1
export annot=$2
export cluster=$3

#######################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

########################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathGenome=${path}/data/genome_sequences/${sp}
export pathScripts=${path}/scripts/gene_families

export release=103
export genomerelease=103

####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${release}
fi

if [ ${annot} = "EnsemblStringTie" ]; then
    export pathAnnot=${pathStringTie}
    export suffix=combined_annotations_StringTie_Ensembl
fi

##############################################################

perl ${pathScripts}/extract.cDNA.sequences.exon.blocks.pl --pathExonBlocks=${pathAnnot}/ExonBlocks_${suffix}.txt  --pathGenomeSequence=${pathGenome}/genome_ensembl${genomerelease}.fa --pathOutput=${pathAnnot}/ExonBlocks_${suffix}_cDNASequences.fa

###############################################################
