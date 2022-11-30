#!/bin/bash

export sp=$1
export cluster=$2

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi


if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTieAssembly=${path}/results/stringtie_assembly/${sp}
export pathScripts=${path}/scripts/transcript_assembly

export release=103

#############################################################################

if [ -e ${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_withTxInfo.gtf ]; then
    echo "file with tx info already there"
else
    perl ${pathScripts}/format.gtf.pl --pathInputGTF=${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl.gtf --pathOutputGTF=${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_withTxInfo.gtf
fi

#############################################################################

cp ${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_withTxInfo.gtf ${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_withTxInfo_withCDSInfo.gtf

grep $'\t'CDS$'\t' ${pathEnsembl}/AllTranscripts_Ensembl103.gtf >> ${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_withTxInfo_withCDSInfo.gtf

#############################################################################
