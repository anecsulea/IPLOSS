#!/bin/bash

export sp=$1
export sample=$2
export annot=$3
export nthreads=$4
export cluster=$5

####################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
    export pathTools=/beegfs/home/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools
fi

####################################################################################

export pathDocs=${path}/docs
export pathRNASeq=${path}/data/RNASeq/${sp}
export pathResults=${path}/results/gene_expression_estimation/${sp}/${annot}
export pathIndexes=${path}/results/kallisto_indexes/${sp}
export pathScripts=${path}/scripts/gene_expression_estimation

####################################################################################

## TruSeq stranded

export library=fr-firststrand
export ss="--rf-stranded"

echo ${sample} ${library} ${ss}

####################################################################################

export file=`grep ${sp} ${pathDocs}/Samples.txt | grep $'\t'${sample}$'\t' | cut -f 7`

if [ -e ${pathRNASeq}/${file} ]; then

    if [ -e ${pathResults}/${sample} ]; then
	echo "path exists"
    else
	mkdir -p ${pathResults}/${sample}
    fi

    if [ -e ${pathResults}/${sample}/abundance.tsv ]; then
	echo "already done"
    else
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}

	if [ ${cluster} = "pbil" ]; then
	    echo "#SBATCH --job-name=kallisto_${sample}" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    echo "#SBATCH --output=${pathScripts}/std_out_kallisto_${sp}_${sample}" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    echo "#SBATCH --error=${pathScripts}/std_err_kallisto_${sp}_${sample}" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    echo "#SBATCH --time=2:00:00" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    echo "#SBATCH --mem=5G" >>  ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	fi

	echo "singularity exec -B ${path} ${pathTools}/kallisto.sif kallisto quant --single -l 200.0 -s 20 --bias -t ${nthreads} ${ss} -o ${pathResults}/${sample} --index ${pathIndexes}/${annot} ${pathRNASeq}/${file} " >> ${pathScripts}/bsub_script_kallisto_${sp}_${sample}

	if [ ${cluster} = "pbil" ]; then
	    sbatch ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	fi

	if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbillocal" ]; then
	    chmod a+x ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	    ${pathScripts}/bsub_script_kallisto_${sp}_${sample}
	fi
    fi
else
    echo "cannot find "${file}
fi

####################################################################################
