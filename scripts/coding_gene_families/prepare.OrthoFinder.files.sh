#!/bin/bash

export dataset=$1
export cluster=$2

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools/OrthoFinder/tools
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
    export pathTools=/beegfs/home/necsulea/Tools/OrthoFinder/tools
fi

export pathGeMoMa=${path}/results/CDS_annotation
export pathEnsemblProteins=${path}/data/protein_sequences/Ensembl103
export pathResults=${path}/results/gene_families/OrthoFinder_${dataset}
export pathScripts=${path}/scripts/coding_gene_families

##########################################################################

if [ -e ${pathResults} ]; then
    echo "output path already there"
else
    mkdir -p ${pathResults}
fi

##########################################################################

if [ ${dataset} = "Ensembl103" ]; then

    for file in `ls  ${pathEnsemblProteins} `
    do
	export sp=`echo ${file} | cut -f 1 -d '.'`
	
	if [ -e ${pathEnsemblProteins}/primary_transcripts/${file} ]; then
	    echo "primary transcripts already done"
	else
	    python ${pathTools}/primary_transcript.py ${pathEnsemblProteins}/${file}
	fi
	
	ln -s ${pathEnsemblProteins}/primary_transcripts/${file} ${pathResults}/${sp}.fa
    done
fi

##########################################################################

if [ ${dataset} = "GeMoMa_pairwise" ]; then
    for sp in Chicken Duck
    do
	export pathAnnot=${pathGeMoMa}/${sp}/GeMoMa/combined
	export file=combined_protein_sequences.faa

	if [ -e ${pathAnnot}/${file} ]; then
	    if [ -e ${pathAnnot}/primary_transcripts/${file} ]; then
		echo "primary transcripts already done"
	    else
		python ${pathTools}/primary_transcript.py ${pathAnnot}/${file}
	    fi
	    
	    ln -s ${pathEnsemblProteins}/primary_transcripts/${file} ${pathResults}/${sp}.fa
	    
	else
	    echo "cannot find protein sequences! "${pathAnnot}/${file} 
	    exit
	fi
    done
fi

##########################################################################

