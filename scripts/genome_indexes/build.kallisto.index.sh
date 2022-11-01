#!/bin/bash

export species=$1
export source=$2
export cluster=$3

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

export pathStringTie=${path}/results/stringtie_assembly/${species}
export pathEnsembl=${path}/data/ensembl_annotations/${species}
export pathResults=${path}/results/kallisto_indexes/${species}

export release=103

####################################################################################

if [ ${source} = "Ensembl" ]; then
    export pathFasta=${pathEnsembl}/AllTranscripts_Ensembl${release}_noMT_norRNA.fa
    export prefix=AllTranscripts_Ensembl${release}
fi

if [ ${source} = "EnsemblStringTie" ]; then
    export pathFasta=${pathStringTie}/combined_annotations_StringTie_Ensembl_noMT_norRNA.fa
    export prefix=EnsemblStringTie
fi

####################################################################################

if [ -e ${pathResults} ]; then
    echo "path results already there"
else
    mkdir -p ${pathResults}
fi

####################################################################################

singularity exec -B ${path} ${pathTools}/kallisto.sif kallisto index -i ${pathResults}/${prefix} ${pathFasta}

####################################################################################
