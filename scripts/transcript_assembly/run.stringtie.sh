#!/bin/bash

export sp=$1
export cluster=$2
export nthreads=$3

#############################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
else
    if [ ${cluster} = "cloud" ]; then
	export path=/ifb/data/mydatalocal/IPLOSS
    else
	echo "unknown cluster"
	exit
    fi
fi

#############################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathAlignments=${path}/results/RNASeq_alignments/${sp}/all_samples
export pathResults=${path}/results/stringtie_assembly/${sp}
export pathScripts=${path}/scripts/transcript_assembly

#############################################################################

## stringtie 2.2.1

#############################################################################

export ensrelease=103

export gtffile=FilteredTranscripts_Ensembl${ensrelease}.gtf

export refpar="-G ${pathEnsembl}/${gtffile}"

#############################################################################

if [ -e ${pathResults} ]; then
    echo "output already there"
else
    mkdir -p ${pathResults}
fi

#############################################################################

if [ -e ${pathResults}/assembled_transcripts.gtf ]&&[ -e ${pathResults}/gene_abundance.txt ]; then
    echo "already done"
else
    if [ -e ${pathAlignments}/accepted_hits.bam ]; then
	echo "#!/bin/bash " > ${pathScripts}/bsub_script_stringtie

	if [ ${cluster} = "pbil" ]; then
	    echo "#SBATCH --job-name=stringtie" >>  ${pathScripts}/bsub_script_stringtie
	    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_stringtie
	    echo "#SBATCH --output=${pathScripts}/std_out_${sp}" >>  ${pathScripts}/bsub_script_stringtie
	    echo "#SBATCH --error=${pathScripts}/std_err_${sp}" >>  ${pathScripts}/bsub_script_stringtie
	    echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script_stringtie ## 1 CPU
	    echo "#SBATCH --time=2:00:00" >>  ${pathScripts}/bsub_script_stringtie ## 2 hours
	    echo "#SBATCH --mem=5G" >>  ${pathScripts}/bsub_script_stringtie ## 4g per CPU
	else
	    echo "#!/bin/bash" > ${pathScripts}/bsub_script_stringtie
	fi

    	echo "stringtie ${pathAlignments}/accepted_hits.bam ${refpar} -m 250 -a 10 -f 0.05 -p ${nthreads} -o ${pathResults}/assembled_transcripts.gtf -A ${pathResults}/gene_abundance.txt --rf -c 5 -s 10 -M 0.5 -g 100 ">> ${pathScripts}/bsub_script_stringtie

	if [ ${cluster} = "pbil" ]; then
	    sbatch ${pathScripts}/bsub_script_stringtie
	fi

	if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbillocal" ]; then
	    chmod a+x ${pathScripts}/bsub_script_stringtie
	    ${pathScripts}/bsub_script_stringtie >& ${pathResults}/std_out_error_stringtie.txt
	fi
    fi
fi
#############################################################################
