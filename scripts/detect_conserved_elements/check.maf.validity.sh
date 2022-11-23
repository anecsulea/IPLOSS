#!/bin/bash

export dataset=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
fi

export pathAln=${path}/results/whole_genome_alignments/${dataset}
export pathScripts=${path}/scripts/detect_conserved_elements

#########################################################################

for file in `ls ${pathAln} | grep maf`
do
    echo ${file} 
    perl ${pathScripts}/check.maf.validity.pl --pathMAF=${pathaln}/${file} 
done

#########################################################################
