#!/bin/bash

export sp=$1
export ref=$2
export sample=$3
export cluster=$4

###################################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export pathStringTie=${path}/results/stringtie_assembly/${sp}/reference_${ref}_${sample}
export pathAlignments=${path}/results/RNASeq_alignments/${sp}/${sample}
export pathScripts=${path}/scripts/transcript_assembly

###################################################################################

export pathGTF=${pathStringTie}/assembled_transcripts.gtf
export pathCorrectJunctions=${pathAlignments}/junctions.txt
export pathWrongJunctions=${pathAlignments}/junctions_wrongstrand.txt

###################################################################################

echo "#!/bin/bash" >  ${pathScripts}/bsub_script_splice

echo "perl ${pathScripts}/check.splice.junctions.pl --pathAnnotGTF=${pathGTF} --pathCorrectJunctions=${pathCorrectJunctions} --pathWrongJunctions=${pathWrongJunctions} --pathOutput=${pathStringTie}/SpliceJunctionsStats.txt"  >>  ${pathScripts}/bsub_script_splice

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_splice
    ${pathScripts}/bsub_script_splice
fi

###################################################################################
