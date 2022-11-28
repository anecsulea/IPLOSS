#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
    export pathTools=/beegfs/home/${USER}/IPLOSS
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/IPLOSS
    export pathTools=/sps/biometr/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools
fi

export pathGeMoMa=${path}/results/CDS_annotation/${sp}/GeMoMa/combined
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathScripts=${path}/scripts/CDS_annotation

#########################################################################

if [ -e ${pathGeMoMa}/filtered_predictions.gtf ]; then
    echo "GTF file already there"
else
    gffread ${pathGeMoMa}/filtered_predictions.gff -T -o  ${pathGeMoMa}/filtered_predictions.gtf 
fi

#########################################################################

perl ${pathScripts}/extract.compatible.CDS.pl --pathCDSAnnotation=${pathGeMoMa}/filtered_predictions.gtf --pathTranscriptAnnotation=${pathStringTie}/combined_annotations_StringTie_Ensembl.gtf --pathOutput=${pathGeMoMa}/compatible_CDS.txt

#########################################################################
