#!/bin/bash

export sp=$1
export annot=$2
export sample=$3

export cluster="cloud"

#####################################################################

if [ ${cluster} = "cloud" ]; then
    export path="/ifb/data/mydatalocal/IPLOSS"
fi

export pathResults=${path}/results/snRNASeq_analysis/${sp}/${annot}/${sample}

#####################################################################

if [ -e ${pathResults}/cellphonedb_meta.txt ]&&[ -e ${pathResults}/cellphonedb_count.txt ]; then
    echo "ok, found input files"
else
    echo "did not find input files"
    echo ${pathResults}/cellphonedb_meta.txt
    echo ${pathResults}/cellphonedb_count.txt
    exit
fi

#####################################################################

conda activate cpdb

cd ${pathResults}

cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt

conda deactivate

#####################################################################
