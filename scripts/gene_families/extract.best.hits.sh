#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3
export aln=$4
export cluster=$5

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

####################################################################################

export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathResults=${path}/results/liftOver_gene_families
export pathScripts=${path}/scripts/gene_families

export release=103

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}
fi

if [ ${annot} = "EnsemblStringTie" ]; then
    export prefix=combined_annotations_StringTie_Ensembl
    export pathExons=${pathStringTie}
fi

####################################################################################

if [ ${sp1} = ${sp2} ]; then
    echo "cannot project from "${sp1}" to "${sp2}
    exit
fi

#################################################################################

perl ${pathScripts}/extract.best.hits.pl --species1=${sp1} --species2=${sp2} --pathExonBlocks1=${pathExons}/${sp1}/ExonBlocks_${prefix}.txt --pathExonBlocks2=${pathExons}/${sp2}/ExonBlocks_${prefix}.txt --pathUTR1=NA --pathUTR2=NA --pathAlignmentStats=${pathResults}/AlignmentStatistics_ExonBlocks_${aln}_${sp1}_${sp2}_${annot}.txt  --minAlignedFraction=0 --minRatioSecondBest=1.1 --pathBestHits12=${pathResults}/BestHits_${aln}_${sp1}_${sp2}_${annot}.txt --pathBestHits21=${pathResults}/BestHits_${aln}_${sp2}_${sp1}_${annot}.txt --pathReciprocalBestHits=${pathResults}/ReciprocalBestHits_${aln}_${sp1}_${sp2}_${annot}.txt

#################################################################################
