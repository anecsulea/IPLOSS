#!/bin/bash

## adapted from script by Menghan Wang

export sp=$1
export sample=$2
export annot=$3
export cluster=$4
export nthreads=$5

####################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathRawData=${path}/data/single_cell_data/rawData
export pathResults=${path}/results/snATACSeq_analysis/${sp}/${annot}
export pathIndex=${path}/results/snATACSeq_indexes/${sp}/${annot}

####################################################################################

if [ -e ${pathResults} ]; then
    echo "output dir already exists"
else
    mkdir -p ${pathResults}
fi

####################################################################################

export pathFastQs=""
export samples=""
export lensuffix=24 ## length of suffix, eg 24 for _S1_L001_I1_001.fastq.gz

for dir in `grep ${sample} ${pathRawData}/samples.txt | cut -f 2`
do
    export pathFastQs=${pathRawData}/${dir},${pathFastQs}
   
    export filename=`ls ${pathRawData}/${dir} | grep I1 | grep fastq.gz`
    export suffix=${filename: -${lensuffix}} 
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

cellranger-atac count --id=${sample} \
                   --reference=${pathIndex} \
                   --fastqs=${pathFastQs} \
                   --sample=${samples} \
                   --localcores=${nthreads} --localmem=60

###########################################################################################################
