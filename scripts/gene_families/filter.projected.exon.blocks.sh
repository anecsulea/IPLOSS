#!/bin/bash

export ref=$1
export tg=$2
export annot=$3
export cluster=$4

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

####################################################################################

export pathProjections=${path}/results/exon_projections
export pathEnsembl=${path}/data/ensembl_annotations
export pathStringTie=${path}/results/stringtie_assembly
export pathScripts=${path}/scripts/gene_families

export release=103

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}/${ref}
fi

if [ ${annot} = "EnsemblStringTie" ]; then
    export prefix=ExonBlocks_combined_annotations_StringTie_Ensembl
    export pathExons=${pathStringTie}/${ref}
fi

####################################################################################

export minsizeratio=0.5
export maxsizeratio=2

####################################################################################

if [ ${ref} = ${tg} ]; then
    echo "cannot project from "${ref}" to "${tg}
    exit
fi

####################################################################################

perl ${pathScripts}/filter.projected.exon.blocks.pl --pathExonBlocks=${pathExons}/${prefix}.txt --pathProjectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}.bed --collapse=100 --minSizeRatio=${minsizeratio} --maxSizeRatio=${maxsizeratio} --pathOutputFilteredExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step1.txt  --pathOutputRejectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_RejectedProjectedExons_Step1.txt --pathOutputLog=${pathScripts}/logs/log_filter_exon_blocks_from${ref}_to${tg}.txt

####################################################################################
