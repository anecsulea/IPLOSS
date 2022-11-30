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
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

####################################################################################

export pathProjections=${path}/results/exon_projections
export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
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

if [ ${ref} = ${tg} ]; then
    echo "cannot project from "${ref}" to "${tg}
    exit
fi

####################################################################################

perl ${pathScripts}/filter.projected.genes.exon.blocks.pl --pathExonBlocks=${pathExons}/${prefix}.txt --pathProjectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step1.txt --maxIntronSizeRatio=100 --maxAddedIntronSize=1000000 --pathSyntenyPredictions=NA --syntenyRange=1000000 --pathOutputLog=${pathScripts}/logs/log_filter_genes_exon_blocks_from${ref}_to${tg}.txt --pathOutputFilteredExons=${pathProjections}/From${ref}_To${tg}_${prefix}_FilteredProjectedExons_Step2.txt --pathOutputRejectedExons=${pathProjections}/From${ref}_To${tg}_${prefix}_RejectedProjectedExons_Step2.txt 

####################################################################################
