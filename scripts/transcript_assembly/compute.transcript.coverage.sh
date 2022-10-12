#!/bin/bash

export sp=$1
export ref=$2
export sample=$3
export cluster=$4

###################################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export pathAlignments=${path}/results/RNASeq_alignments/${sp}/${sample}
export pathResults=${path}/results/stringtie_assembly/${sp}/reference_${ref}_${sample}
export pathScripts=${path}/scripts/transcript_assembly

###############################################################################

export pathGTF=${pathResults}/assembled_transcripts.gtf

if [ -e ${pathResults}/coverage_sense_antisense ]; then
    echo "dir output exists"
else
    mkdir -p ${pathResults}/coverage_sense_antisense
fi

###############################################################################

if [ -e ${pathResults}/coverage_sense_antisense/CoverageExons.txt ]&&[ -e ${pathResults}/coverage_sense_antisense/CoverageTranscripts.txt ]; then
    echo "already done"
else

    if [ -e ${pathAlignments}/coverage_reverse.bedGraph.gz ]&&[ -e ${pathAlignments}/coverage_forward.bedGraph.gz ]; then
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_coverage

	echo "perl ${pathScripts}/compute.transcript.coverage.pl --pathAnnotGTF=${pathGTF} --pathCoverageForward=${pathAlignments}/coverage_forward.bedGraph.gz --pathCoverageReverse=${pathAlignments}/coverage_reverse.bedGraph.gz --pathOutputExons=${pathResults}/coverage_sense_antisense/CoverageExons.txt --pathOutputTranscripts=${pathResults}/coverage_sense_antisense/CoverageTranscripts.txt"  >>  ${pathScripts}/bsub_script_coverage


	if [ ${cluster} = "cloud" ]; then
	    chmod a+x ${pathScripts}/bsub_script_coverage
	    ${pathScripts}/bsub_script_coverage
	fi

    fi
fi

###################################################################################
