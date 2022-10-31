#!/bin/bash

export sp=$1
export sample=$2
export cluster=$3

###################################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export pathAlignments=${path}/results/RNASeq_alignments/${sp}/${sample}
export pathScripts=${path}/scripts/RNASeq_alignments

###################################################################################

## we have rf libraries, bedtools doesn't manage this well

###################################################################################

if [ -e ${pathAlignments}/coverage_reverse.bedGraph.gz ]&&[ -e ${pathAlignments}/coverage_forward.bedGraph.gz ]; then
    echo "already done"
else
    echo "#!/bin/bash" >  ${pathScripts}/bsub_script_coverage_${sp}

    ## forward strand

    echo "bedtools genomecov -strand \"-\" -bg -split -ibam ${pathAlignments}/accepted_hits.bam > ${pathAlignments}/coverage_forward.bedGraph" >>  ${pathScripts}/bsub_script_coverage_${sp}

    echo "gzip ${pathAlignments}/coverage_forward.bedGraph" >>  ${pathScripts}/bsub_script_coverage_${sp}

    ## reverse strand

    echo "bedtools genomecov -strand \"+\" -bg -split -ibam ${pathAlignments}/accepted_hits.bam > ${pathAlignments}/coverage_reverse.bedGraph" >>  ${pathScripts}/bsub_script_coverage_${sp}

    echo "gzip ${pathAlignments}/coverage_reverse.bedGraph" >>  ${pathScripts}/bsub_script_coverage_${sp}

    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_coverage_${sp}
	${pathScripts}/bsub_script_coverage_${sp}
    fi
fi

###################################################################################
