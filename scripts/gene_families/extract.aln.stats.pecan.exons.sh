#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3
export cluster=$4

####################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathResults=${path}/results/liftOver_gene_families
export pathScripts=${path}/scripts/gene_families

export release=103

####################################################################

if [ ${annot} = "EnsemblStringTie" ]; then
    export prefixExons=ExonBlocks_combined_annotations_StringTie_Ensembl
fi

if [ ${annot} = "Ensembl" ]; then
    export prefixExons=ExonBlocks_FilteredTranscripts_Ensembl${release}
fi

####################################################################
        
perl ${pathScripts}/extract.aln.stats.pecan.exons.pl --species1=${sp1} --species2=${sp2} --pathClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_${prefixExons}.txt --dirPecan=${pathResults}/pecan_alignments_projection_clusters/${sp1}_${sp2}_${annot}/ --minAlignmentLength=0 --pathOutput=${pathResults}/AlignmentStatistics_ExonBlocks_Pecan_${sp1}_${sp2}_${annot}.txt 

####################################################################
