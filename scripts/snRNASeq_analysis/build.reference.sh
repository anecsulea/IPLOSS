#!/bin/bash

## adapted from script by Menghan Wang

export sp=$1
export annot=$2
export cluster=$3
export nthreads=$4

####################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
    export pathTools=/beegfs/home/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools
fi

####################################################################################

export pathGenome=${path}/data/genome_sequences/${sp}
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathIndexes=${path}/results/snRNA_indexes/${sp}

export release=103

## cell ranger version 7.0.1

####################################################################################

if [ -e ${pathIndexes} ]; then
    echo "path output already there"
else
    mkdir -p ${pathIndexes} 
fi

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf
fi

####################################################################################

if [ ${annot} = "EnsemblStringTie" ]; then
    export pathGTF=${pathStringTie}/combined_annotations_StringTie_Ensembl.gtf
fi

####################################################################################

cd ${pathIndexes}

cellranger mkref --genome=${annot}  --fasta=${pathGenome}/genome_ensembl${release}.fa --genes=${pathGTF} --nthreads=${nthreads}

########################################################################
