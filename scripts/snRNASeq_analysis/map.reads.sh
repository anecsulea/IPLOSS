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

## command each sample
#sbatch step1_rna_mappingPipe.sh -i chicken_wholeGT_male_rna -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/chicken_gg6/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603041140936-60870697,/scicore/projects/openbis/userstore/duw_tschopp/20220603041230793-60870699 -s BSSE_QGF_206707_HCCHMDRX2_1_Gg_W_GT_male_33_scRNA_SI-TT-H1_,BSSE_QGF_206707_HCCHMDRX2_2_Gg_W_GT_male_33_scRNA_SI-TT-H1_
#sbatch step1_rna_mappingPipe.sh -i chicken_wholeGT_male_rna_ensembl -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/chicken_gg6_ensembl/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603041140936-60870697,/scicore/projects/openbis/userstore/duw_tschopp/20220603041230793-60870699 -s BSSE_QGF_206707_HCCHMDRX2_1_Gg_W_GT_male_33_scRNA_SI-TT-H1_,BSSE_QGF_206707_HCCHMDRX2_2_Gg_W_GT_male_33_scRNA_SI-TT-H1_
#sbatch step1_rna_mappingPipe_withintron.sh -i chicken_wholeGT_male_rna_withintron -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/chicken_gg6/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603041140936-60870697,/scicore/projects/openbis/userstore/duw_tschopp/20220603041230793-60870699 -s BSSE_QGF_206707_HCCHMDRX2_1_Gg_W_GT_male_33_scRNA_SI-TT-H1_,BSSE_QGF_206707_HCCHMDRX2_2_Gg_W_GT_male_33_scRNA_SI-TT-H1_

#sbatch step1_rna_mappingPipe.sh -i duck_wholeGT_male_rna -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/duck_ap1/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603041254383-60870700,/scicore/projects/openbis/userstore/duw_tschopp/20220603041203049-60870698 -s BSSE_QGF_206706_HCCHMDRX2_2_Ap_W_GT_male_33_scRNA_SI-TT-G12_,BSSE_QGF_206706_HCCHMDRX2_1_Ap_W_GT_male_33_scRNA_SI-TT-G12_
#sbatch step1_rna_mappingPipe.sh -i duck_wholeGT_male_rna_ensembl -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/duck_ap1_ensembl/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603041254383-60870700,/scicore/projects/openbis/userstore/duw_tschopp/20220603041203049-60870698 -s BSSE_QGF_206706_HCCHMDRX2_2_Ap_W_GT_male_33_scRNA_SI-TT-G12_,BSSE_QGF_206706_HCCHMDRX2_1_Ap_W_GT_male_33_scRNA_SI-TT-G12_
#sbatch step1_rna_mappingPipe_withintron.sh -i duck_wholeGT_male_rna_withintron -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/duck_ap1/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603041254383-60870700,/scicore/projects/openbis/userstore/duw_tschopp/20220603041203049-60870698 -s BSSE_QGF_206706_HCCHMDRX2_2_Ap_W_GT_male_33_scRNA_SI-TT-G12_,BSSE_QGF_206706_HCCHMDRX2_1_Ap_W_GT_male_33_scRNA_SI-TT-G12_

#sbatch step1_rna_mappingPipe.sh -i duck_wholeGT_female_rna_ensembl -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/duck_ap1_ensembl/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603040844650-60870686,/scicore/projects/openbis/userstore/duw_tschopp/BSSE_QGF_206697_HCCHMDRX2_2 -s BSSE_QGF_206697_HCCHMDRX2_1_Ap_W_GT_female_33_scRNA_SI-TT-F10_,BSSE_QGF_206697_HCCHMDRX2_2_Ap_W_GT_female_33_scRNA_SI-TT-F10_
#sbatch step1_rna_mappingPipe.sh -i duck_wholeGT_female_rna -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/duck_ap1/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603040844650-60870686,/scicore/projects/openbis/userstore/duw_tschopp/BSSE_QGF_206697_HCCHMDRX2_2 -s BSSE_QGF_206697_HCCHMDRX2_1_Ap_W_GT_female_33_scRNA_SI-TT-F10_,BSSE_QGF_206697_HCCHMDRX2_2_Ap_W_GT_female_33_scRNA_SI-TT-F10_
#sbatch step1_rna_mappingPipe_withintron.sh -i duck_wholeGT_female_rna_withintron -t /scicore/home/tschoppp/GROUP/references/genomes/forGenitalia/duck_ap1/ -f /scicore/projects/openbis/userstore/duw_tschopp/20220603040844650-60870686,/scicore/projects/openbis/userstore/duw_tschopp/BSSE_QGF_206697_HCCHMDRX2_2 -s BSSE_QGF_206697_HCCHMDRX2_1_Ap_W_GT_female_33_scRNA_SI-TT-F10_,BSSE_QGF_206697_HCCHMDRX2_2_Ap_W_GT_female_33_scRNA_SI-TT-F10_




