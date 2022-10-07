#!/bin/bash

export sp=$1
export cluster=$2

##################################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export pathAnnot=${path}/data/ensembl_annotations/${sp}
export pathScripts=${path}/scripts/process_ensembl_annotations

export release=103

##################################################################################################

perl ${pathScripts}/extract.readthrough.transcripts.pl --pathExonCoords=${pathAnnot}/ExonCoords_Ensembl${release}.txt  --pathExonsTranscripts=${pathAnnot}/ExonsTranscripts_Ensembl${release}.txt --pathTranscriptInfo=${pathAnnot}/TranscriptInfo_Ensembl${release}.txt --pathGeneInfo=${pathAnnot}/GeneInfo_Ensembl${release}.txt --monoexonicBiotypes="Mt_tRNA,Mt_rRNA,IG_D_gene,IG_J_gene,snoRNA,misc_RNA,miRNA,snRNA,rRNA"  --pathOutput=${pathAnnot}/ReadthroughTranscripts_Ensembl${release}.txt

##################################################################################################
