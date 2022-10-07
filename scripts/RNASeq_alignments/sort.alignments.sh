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

export pathResults=${path}/results/RNASeq_alignments/${sp}
export pathScripts=${path}/scripts/RNASeq_alignments

#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.bam ]; then
    echo "already done (or ongoing)"
else
    if [ -e ${pathResults}/${sample}/accepted_hits.sam ]; then
	echo "sorting alignments"

	echo "#!/bin/bash" > ${pathScripts}/bsub_script_sort_${sp}_${sample}
	
	if [ ${cluster} = "pbil" ]; then
	    echo "#!/bin/bash " > ${pathScripts}/bsub_script_sort_${sp}_${sample}
	    echo "#SBATCH --job-name=sort_${sample}" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	    echo "#SBATCH --output=${pathScripts}/std_out_${sample}" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	    echo "#SBATCH --error=${pathScripts}/std_err_${sample}" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	    echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample} ## ${nthreads} CPU
	    echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample} ## 24 hours
	    echo "#SBATCH --mem=10G" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample} ## 5g per CPU
	fi	
	
	echo "samtools sort -m 10G -@ ${nthreads} -O bam -o ${pathResults}/${sample}/accepted_hits.bam ${pathResults}/${sample}/accepted_hits.sam" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample} 
	
	echo "if [ -e ${pathResults}/${sample}/accepted_hits.bam ]; then " >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	echo "rm ${pathResults}/${sample}/accepted_hits.sam" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	echo "else" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	echo "echo 'not done'" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	echo "fi" >>  ${pathScripts}/bsub_script_sort_${sp}_${sample}
	
	if [ ${cluster} = "cloud" ]; then
	    chmod a+x ${pathScripts}/bsub_script_sort_${sp}_${sample}
	    ${pathScripts}/bsub_script_sort_${sp}_${sample}
	fi
	
	if [ ${cluster} = "pbil" ]; then
	    sbatch ${pathScripts}/bsub_script_sort_${sp}_${sample}
	fi
    else
	echo "cannot find sam file"
    fi
fi

#############################################################################
