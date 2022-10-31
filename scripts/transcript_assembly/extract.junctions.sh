#!/bin/bash

export sp=$1
export sample=$2
export cluster=$3

######################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export ensrelease=103

export pathAlignments=${path}/results/RNASeq_alignments/${sp}/${sample}
export pathSequence=${path}/data/genome_indexes/${sp}/genome_ensembl${ensrelease}.fa
export pathScripts=${path}/scripts/transcript_assembly

######################################################################################

if [ -e ${pathAlignments}/junctions.txt ]; then
    echo "already done"
else

    if [ -e ${pathAlignments}/accepted_hits.bam ]; then
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_junctions_${sp}

	if [ ${cluster} = "pbil" ]; then
	    echo "#!/bin/bash " > ${pathScripts}/bsub_script_junctions_${sp}
	    echo "#SBATCH --job-name=junctions_${sp}" >>  ${pathScripts}/bsub_script_junctions_${sp}
	    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_junctions_${sp}
	    echo "#SBATCH --output=${pathScripts}/std_out_${sp}" >>  ${pathScripts}/bsub_script_junctions_${sp}
	    echo "#SBATCH --error=${pathScripts}/std_err_${sp}" >>  ${pathScripts}/bsub_script_junctions_${sp}
	    echo "#SBATCH --cpus-per-task=1" >>  ${pathScripts}/bsub_script_junctions_${sp} ## 1 CPU
	    echo "#SBATCH --time=1:00:00" >>  ${pathScripts}/bsub_script_junctions_${sp} ## 4 hours
	    echo "#SBATCH --mem=8G" >>  ${pathScripts}/bsub_script_junctions_${sp} ## 8g per CPU
	fi

	echo "perl ${pathScripts}/extract.junctions.pl --pathAln=${pathAlignments}/accepted_hits.bam --pathGenomeSequence=${pathSequence} --anchordetection=8 --anchorquantification=5 --maxmismatch=0.02  --pathOutput=${pathAlignments}/junctions.txt --pathOutputWrongStrand=${pathAlignments}/junctions_wrongstrand.txt " >> ${pathScripts}/bsub_script_junctions_${sp}

	if [ ${cluster} = "cloud" ]||[ ${cluster} = "pbillocal" ]; then
	    chmod a+x ${pathScripts}/bsub_script_junctions_${sp}
	    ${pathScripts}/bsub_script_junctions_${sp}
	fi

	if [ ${cluster} = "pbil" ]; then
	    sbatch ${pathScripts}/bsub_script_junctions_${sp}
	fi
    else
	echo "alignment is not there"
    fi
fi


######################################################################################
######################################################################################
