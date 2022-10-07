#!/bin/bash

export sp=$1
export sample=$2
export cluster=$3
export nthreads=$4

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

export pathDocs=${path}/docs
export pathRNASeq=${path}/data/RNASeq/${sp}
export pathIndexes=${path}/data/genome_indexes/${sp}
export pathResults=${path}/results/RNASeq_alignments/${sp}
export pathScripts=${path}/scripts/RNASeq_alignments

## hisat2-2.2.1

#############################################################################
 
export prefix=genome_ensembl103

export pathIndex=${pathIndexes}/${prefix}

#############################################################################

if [ -e ${pathResults}/${sample} ]; then
    echo "output dir exists"
else
    mkdir -p ${pathResults}/${sample}
fi

#############################################################################

if [ ${sp} = "Chicken" ]||[ ${sp} = "Duck" ]; then
    export strand="--rna-strandness R" ## TruSeq
else
    export strand=""
fi

#############################################################################

export fastqfile=`grep ${sp} ${pathDocs}/Samples.txt | grep $'\t'${sample}$'\t' | cut -f 7`
export seqtype="single_end"

echo ${sample} ${fastqfile}

if [ -e ${pathRNASeq}/${fastqfile} ]; then
    echo "found file for "${sp} ${sample}
else
    echo "cannot find data for "${sp} ${sample}
    exit
fi
  
#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.sam.gz ]||[ -e ${pathResults}/${sample}/accepted_hits.bam ]||[ -e ${pathResults}/${sample}/accepted_hits.sam ]||[ -e ${pathResults}/${sample}/accepted_hits_clean.sam ]||[ -e ${pathResults}/${sample}/accepted_hits_clean.sam.gz ]; then
    echo "already done (or ongoing)"
else
    
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_hisat_${sp}_${sample}

    if [ ${cluster} = "pbil" ]; then
	echo "#!/bin/bash " > ${pathScripts}/bsub_script_hisat_${sp}_${sample}
	echo "#SBATCH --job-name=hisat_${sample}" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample}
	echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample}
	echo "#SBATCH --output=${pathScripts}/std_out_${sample}" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample}
	echo "#SBATCH --error=${pathScripts}/std_err_${sample}" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample}
	echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample} ## ${nthreads} CPU
	echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample} ## 24 hours
	echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_hisat_${sp}_${sample} ## 5g per CPU
    fi	

    if [ ${seqtype} = "single_end" ]; then
	echo "single-end"
	echo "hisat2 --seed 19 -p ${nthreads} -x ${pathIndex} -U ${pathRNASeq}/${fastqfile} -S ${pathResults}/${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/${sample}/novel_splicesites.txt >& ${pathResults}/${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat_${sp}_${sample}
    fi
        
    if [ ${seqtype} = "paired_end" ]; then
	echo "not sure what to do"
	exit
    fi
    
    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_hisat_${sp}_${sample}
	${pathScripts}/bsub_script_hisat_${sp}_${sample}
    fi
    
    if [ ${cluster} = "pbil" ]; then
	sbatch ${pathScripts}/bsub_script_hisat_${sp}_${sample}
    fi
    
fi

#############################################################################
