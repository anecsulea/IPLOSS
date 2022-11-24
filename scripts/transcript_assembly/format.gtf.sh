#!/bin/bash

export sp=$1
export cluster=$2

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi


if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTieAssembly=${path}/results/stringtie_assembly/${sp}
export pathScripts=${path}/scripts/transcript_assembly

#############################################################################

perl ${pathScripts}/format.gtf.pl --pathInputGTF=${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl.gtf --pathOutputGTF=${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_withTxInfo.gtf

#############################################################################
