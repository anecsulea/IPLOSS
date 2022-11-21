#!/bin/bash

export ref=$1
export source=$2
export cluster=$3

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

#########################################################################

if [ ${source} = "Reptiles_Ensembl103" ]; then
    export pathAnnotations=${path}/data/ensembl_annotations/${source}
fi

if [ ${source} = "NCBI" ]; then
    export pathAnnotations=${path}/data/NCBI_annotations
fi

export pathScripts=${path}/scripts/CDS_annotation

#########################################################################

if [ -e ${pathAnnotations}/parts ]; then
    echo "dir output exists"
else
    mkdir -p ${pathAnnotations}/parts
fi

#########################################################################

export annotfile=`ls ${pathAnnotations} | grep ${ref}'\.' | grep gff`

echo "annotation file ".${annotfile}." prefix "${ref}

#########################################################################

perl ${pathScripts}/split.GFF.pl --pathGFF=${pathAnnotations}/${annotfile} --maxChrSize=50000000 --dirOutput=${pathAnnotations}/parts --prefixOutput=${ref}

#########################################################################
