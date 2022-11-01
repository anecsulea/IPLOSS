#!/bin/bash

export sp=$1
export annot=$2
export cluster=$3

##############################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
    export pathTools=/beegfs/home/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
    export pathTools=/home/ubuntu/data/mydatalocal/Tools
fi

###########################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathGenome=${path}/data/genome_sequences/${sp}
export pathScripts=${path}/scripts/genome_indexes

export release=103

##############################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}
    export suffix=AllTranscripts_Ensembl${release}
fi

##############################################################

if [ ${annot} = "EnsemblStringTie" ]; then
    export pathGTF=${pathStringTie}
    export suffix=combined_annotations_StringTie_Ensembl
fi

##############################################################

perl ${pathScripts}/extract.cDNA.sequences.pl --pathAnnotGTF=${pathGTF}/${suffix}.gtf --forbiddenChromo=NA --pathGeneInfo=${pathEnsembl}/GeneInfo_Ensembl${release}.txt --forbiddenBiotypes=Mt_rRNA,rRNA,rRNA_pseudogene --pathGenomeSequence=${pathGenome}/genome_ensembl${release}.fa --pathOutput=${pathGTF}/${suffix}_norRNA.fa
   
###############################################################
