#!/bin/bash

## adapted from script by Menghan Wang

export sp=$1
export sample=$2
export lensuffix=$3 
export annot=$4
export cluster=$5
export nthreads=$6

####################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathRawData=${path}/data/single_cell_data/rawData
export pathResults=${path}/results/snRNASeq_analysis/${sp}/${annot}
export pathIndex=${path}/results/snRNASeq_indexes/${sp}/${annot}

## cell ranger version 7.0.1

####################################################################################

if [ -e ${pathResults} ]; then
    echo "output dir already exists"
else
    mkdir -p ${pathResults}
fi

####################################################################################

export pathFastQs=""
export samples=""

for dir in `grep ${sample} ${pathRawData}/samples.txt | cut -f 2`
do
    export pathFastQs=${pathRawData}/${dir},${pathFastQs}
   
    export filename=`ls ${pathRawData}/${dir} | grep I1 | grep fastq.gz`
    export suffix=${filename: -${lensuffix}} ## length of suffix, eg 24 for _S1_L001_I1_001.fastq.gz, 25 otherwise
    export prefix=`basename ${filename} ${suffix}`
    echo ${prefix} ${suffix}
    export samples=${prefix},${samples}
done

export pathFastQs=${pathFastQs::-1}
export samples=${samples::-1}

echo ${pathFastQs}
echo ${samples}

###########################################################################################################

cd ${pathResults}

cellranger count --id=${sample} \
           --transcriptome=${pathIndex} \
           --fastqs=${pathFastQs} \
           --sample=${samples} \
           --include-introns true \
	   --localcores=${nthreads} --localmem=32

###########################################################################################################
