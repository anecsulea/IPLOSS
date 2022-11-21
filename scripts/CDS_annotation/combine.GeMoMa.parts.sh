#!/bin/bash

export target=$1
export ref=$2
export source=$3
export cluster=$4

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export release=103

export pathAllGenomes=${path}/data/genome_sequences
export pathGenomeSequence=${path}/data/genome_sequences/${target}/


if [ ${source} = "Reptiles_Ensembl103_parts" ]; then
    export pathSourceAnnotations=${path}/data/ensembl_annotations/Reptiles_Ensembl103/parts
    export pathSourceGenomes=${path}/data/genome_sequences/Reptiles_Ensembl103
fi

if [ ${source} = "NCBI_parts" ]; then
    export pathSourceAnnotations=${path}/data/NCBI_annotations/parts
    export pathSourceGenomes=${path}/data/genome_sequences/NCBI
fi


export pathResults=${path}/results/CDS_annotation/${target}/GeMoMa
export pathScripts=${path}/scripts/CDS_annotation

#########################################################################

if [ -e ${pathResults}/${ref} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}/${ref}
fi

#########################################################################

if [ -e ${pathResults}/${ref}/final_annotation.gff ]; then
    echo "combined annotations exist, not doing anything"
    exit
fi

#########################################################################

for part in `ls ${pathSourceAnnotations} | grep ${ref} | cut -f 2 -d '.'`
do
    if [ -e ${pathResults}/${ref}_${part}/final_annotation.gff ]; then
	echo ${part}
	
	if [ -e ${pathResults}/${ref}/final_annotation.gff ]; then
	    cat ${pathResults}/${ref}_${part}/final_annotation.gff | grep -v SOFTWARE | grep -v gff-version >> ${pathResults}/${ref}/final_annotation.gff
	else
	    cp ${pathResults}/${ref}_${part}/final_annotation.gff ${pathResults}/${ref}/final_annotation.gff
	fi
    else
	export nbcds=`grep -c CDS ${pathSourceAnnotations}/${ref}.${part}.gff`

	if [ ${nbcds} == '0' ]; then
	    echo "cannot find results for "${part}", but there are no CDS"
	else
	    echo "cannot find results for "${part}", there are "${nbcds}" CDS"
	    exit
	fi
    fi
done

#########################################################################
