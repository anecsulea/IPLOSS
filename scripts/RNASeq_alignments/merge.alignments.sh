#!/bin/bash

export sp=$1
export cluster=$2
export nthreads=$3

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

export pathResults=${path}/results/RNASeq_alignments/${sp}
export pathScripts=${path}/scripts/RNASeq_alignments

#############################################################################

if [ -e ${pathResults}/all_samples ]; then
    echo "output dir exists"
else
    mkdir -p ${pathResults}/all_samples
fi

#############################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_merge_${sp}

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=merge_${sp}" >> ${pathScripts}/bsub_script_merge_${sp}
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_merge_${sp}
    echo "#SBATCH --output=${pathScripts}/std_out_${sample}" >> ${pathScripts}/bsub_script_merge_${sp}
    echo "#SBATCH --error=${pathScripts}/std_err_${sample}" >> ${pathScripts}/bsub_script_merge_${sp}
    echo "#SBATCH --cpus-per-task=${nthreads}" >> ${pathScripts}/bsub_script_merge_${sp}
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_merge_${sp}
    echo "#SBATCH --mem=4G" >> ${pathScripts}/bsub_script_merge_${sp}
fi

echo -n "samtools merge --threads ${nthreads} -o ${pathResults}/all_samples/accepted_hits.bam " >> ${pathScripts}/bsub_script_merge_${sp}

for sample in `ls ${pathResults} | grep -v samples `
do
    if [ -e ${pathResults}/${sample}/accepted_hits.bam ]; then
	echo -n "${pathResults}/${sample}/accepted_hits.bam " >> ${pathScripts}/bsub_script_merge_${sp}
    fi
done

echo "" >> ${pathScripts}/bsub_script_merge_${sp}

#############################################################################

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_merge_${sp}
    ${pathScripts}/bsub_script_merge_${sp}
fi

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    sbatch ${pathScripts}/bsub_script_merge_${sp}
fi

#############################################################################
