#!/bin/bash

export sp=$1
export annot=$2
export cluster=$3

#####################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathScripts=${path}/scripts/process_annotations

export ensrelease=103

#####################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathAnnot=${pathEnsembl}
    export suffix=FilteredTranscripts_Ensembl${ensrelease}
fi

if [ ${annot} = "EnsemblStringTie" ]; then
    export pathAnnot=${pathStringTie}
    export suffix=combined_annotations_StringTie_Ensembl
fi

#####################################################################

perl ${pathScripts}/make.exon.blocks.gtf.pl --pathGTF=${pathAnnot}/${suffix}.gtf --collapseDistance=0  --pathOutputExonBlocks=${pathAnnot}/ExonBlocks_${suffix}.txt

#################################################################################
