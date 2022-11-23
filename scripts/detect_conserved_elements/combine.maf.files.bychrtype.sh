#!/bin/bash

## original script by Alexandre Laverré

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

export alnprefix="aln"

#########################################################################

export firstfile="${alnprefix}_1.maf"

cp ${pathAln}/${firstfile} ${pathAln}/${alnprefix}_macro_chromosomes.maf

for chr in {2..8}
do
    if [ -e ${pathAln}/${alnprefix}_${chr}.maf ]; then
	grep -v "#" ${pathAln}/${alnprefix}_${chr}.maf | sed '1d' >> ${pathHAL}/${alnprefix}_macro_chromosomes.maf
    else
	echo "cannot find file "${pathAln}/${alnprefix}_${chr}.maf 
    fi
done

#########################################################################

export firstfile="${alnprefix}_chr9.maf"

cp ${pathAln}/${firstfile} ${pathHAL}/${alnprefix}_micro_chromosomes.maf

for chr in {10..29}
do
    if [ -e ${pathAln}/${alnprefix}_${chr}.maf ]; then
	grep -v "#" ${pathAln}/${alnprefix}_${chr}.maf | sed '1d' >> ${pathHAL}/${alnprefix}_micro_chromosomes.maf
    else
	echo "cannot find file "${pathAln}/${alnprefix}_${chr}.maf 
    fi
done

#########################################################################

export firstfile="${alnprefix}_chrZ.maf"

cp ${pathAln}/${firstfile} ${pathHAL}/${alnprefix}_sex_chromosomes.maf

## no W chromosome

#########################################################################

export firstfile="${alnprefix}_small.maf"

cp ${pathAln}/${firstfile} ${pathHAL}/${alnprefix}_scaffolds.maf

for file in `ls ${pathAln} | grep PED | grep -v ${firstfile}`
do
    if [ -e ${pathAln}/${file} ]; then
	grep -v "#" ${pathAln}/${file} | sed '1d' >> ${pathHAL}/${alnprefix}_scaffolds.maf
    else
	echo "cannot find file "${pathAln}/${file}
    fi
done

#########################################################################
