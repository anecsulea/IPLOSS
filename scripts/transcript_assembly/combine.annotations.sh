#!/bin/bash

export sp=$1
export ref=$2
export sample=$3
export cluster=$4

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/IPLOSS
fi

export pathEnsembl=${path}/data/ensembl_annotations/${sp}
export pathStringTieAssembly=${path}/results/stringtie_assembly/${sp}/reference_${ref}_${sample}
export pathScripts=${path}/scripts/transcript_assembly

#############################################################################

export ensrelease=103

export ensprefix=FilteredTranscripts_Ensembl${ensrelease} ## no read-through, no monoexonic - we add them later
export stringtieprefix=assembled_transcripts

#############################################################################

perl ${pathScripts}/select.new.transcripts.pl --pathAnnot1=${pathEnsembl}/${ensprefix}.gtf --pathAnnot2=${pathStringTieAssembly}/${stringtieprefix}.gtf --chrList=NA --pathSenseAntisenseCoverage=${pathStringTieAssembly}/coverage_sense_antisense/CoverageTranscripts.txt --minSenseAntisenseRatio=0.05 --pathSpliceJunctionStatistics=${pathStringTieAssembly}/SpliceJunctionsStats.txt --pathOutputDecision=${pathStringTieAssembly}/selected_transcripts_vs_Ensembl.txt

#############################################################################

perl ${pathScripts}/combine.annotations.pl --pathAnnot1=${pathEnsembl}/${ensprefix}.gtf --pathAnnot2=${pathStringTieAssembly}/${stringtieprefix}.gtf --chrList=NA --pathSelectedTranscripts=${pathStringTieAssembly}/selected_transcripts_vs_Ensembl.txt --pathOutputGTF=${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_without_monoexonic.gtf

#############################################################################

cat ${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl_without_monoexonic.gtf ${pathEnsembl}/monoexonicTranscripts_Ensembl${ensrelease}.gtf > ${pathStringTieAssembly}/combined_annotations_StringTie_Ensembl.gtf

#############################################################################
