#!/bin/bash

## adapted from script by Menghan Wang

export sp=$1
export annot=$2
export cluster=$3

####################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
    export pathTools=/beegfs/home/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools
fi

export pathTFMotifs=${path}/data/TF_motifs
export pathGenome=${path}/data/genome_sequences/${sp}
export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathIndexes=${path}/results/snATACSeq_indexes/${sp}

export release=103

## cellranger-atac version 2.1.0 

####################################################################################

if [ -e ${pathIndexes} ]; then
    echo "path output already there"
else
    mkdir -p ${pathIndexes} 
fi

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export pathGTF=${pathEnsembl}/AllTranscripts_Ensembl${release}.gtf
fi

####################################################################################

if [ ${annot} = "EnsemblStringTie" ]; then
    export pathGTF=${pathStringTie}/combined_annotations_StringTie_Ensembl_withTxInfo.gtf
fi

####################################################################################

echo "{" > ${pathIndexes}/config_${annot}.txt
echo "     organism: \""${sp}"\"" >> ${pathIndexes}/config_${annot}.txt
echo "     genome: [\""${annot}"\"]"  >> ${pathIndexes}/config_${annot}.txt
echo "     input_fasta: [\""${pathGenome}/genome_ensembl${release}.fa"\"]"  >> ${pathIndexes}/config_${annot}.txt
echo "     input_gtf: [\""${pathGTF}"\"]"  >> ${pathIndexes}/config_${annot}.txt
if [ ${sp} = "Chicken" ]; then
    echo "     non_nuclear_contigs: [\"MT\"]" >> ${pathIndexes}/config_${annot}.txt
fi
echo "     input_motifs: \""${pathTFMotifs}/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"\"" >> ${pathIndexes}/config_${annot}.txt
echo "}" >> ${pathIndexes}/config_${annot}.txt

####################################################################################

cd ${pathIndexes}

cellranger-atac mkref --config=config_${annot}.txt

####################################################################################
