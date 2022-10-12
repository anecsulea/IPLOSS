#!/bin/bash

export sp=$1
export mode=$2
export cluster=$3

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTieAssembly=${path}/results/stringtie_assembly/${sp}/${mode}
export pathScripts=${path}/scripts/transcript_assembly

#############################################################################

export ensrelease=103

export ensprefix=AllTranscripts_Ensembl${ensrelease} 
export stringtieprefix=assembled_transcripts

#############################################################################

perl ${pathScripts}/compare.annotations.pl --pathAnnot1=${pathStringTieAssembly}/${stringtieprefix}.gtf --pathAnnot2=${pathEnsembl}/${ensprefix}.gtf --pathOutput=${pathStringTieAssembly}/assembled_transcripts_vs_Ensembl.txt

perl ${pathScripts}/compare.annotations.pl --pathAnnot1=${pathEnsembl}/${ensprefix}.gtf --pathAnnot2=${pathStringTieAssembly}/${stringtieprefix}.gtf --pathOutput=${pathStringTieAssembly}/Ensembl_vs_assembled_transcripts.txt

#############################################################################
