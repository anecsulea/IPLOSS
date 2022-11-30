#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3
export cluster=$4

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

####################################################################################

export pathProjections=${path}/results/liftOver_gene_families
export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathResults=${path}/results/liftOver_gene_families
export pathScripts=${path}/scripts/gene_families

export release=103

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}
fi

if [ ${annot} = "EnsemblStringTie" ]; then
    export prefix=ExonBlocks_combined_annotations_StringTie_Ensembl
    export pathExons=${pathStringTie}
fi

####################################################################################

if [ ${sp1} = ${sp2} ]; then
    echo "cannot project from "${sp1}" to "${sp2}
    exit
fi

####################################################################################

perl ${pathScripts}/extract.projection.overlaps.pl --species1=${sp1} --species2=${sp2} --pathExonBlocks1=${pathExons}/${sp1}/${prefix}.txt --pathExonBlocks2=${pathExons}/${sp2}/${prefix}.txt --pathProjectedExons12=${pathProjections}/From${sp1}_To${sp2}_${prefix}_FilteredProjectedExons_Step2.txt --pathProjectedExons21=${pathProjections}/From${sp2}_To${sp1}_${prefix}_FilteredProjectedExons_Step2.txt --pathProjectionMap12=${pathResults}/ProjectionMap_From${sp1}_To${sp2}_${prefix}.txt --pathProjectionMap21=${pathResults}/ProjectionMap_From${sp2}_To${sp1}_${prefix}.txt --pathGeneClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_${prefix}.txt

####################################################################################
