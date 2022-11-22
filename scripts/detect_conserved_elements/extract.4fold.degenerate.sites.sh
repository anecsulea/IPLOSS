#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export refsp=$1
export dataset=$2
export cluster=$3

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathAnnot=${path}/data/ensembl_annotations/${refsp}/by_chr ## annotations for reference species
export pathAlignments=${path}/results/whole_genome_alignments/${dataset} ## alignments in maf format
export pathResults=${path}/results/conserved_elements/${dataset}/4d/

######################################################################

export alnprefix="aln" ## alignment file name
export annotprefix="AllTranscripts_Ensembl103.CDS" ## CDS coordinates for reference species

## extract 4-fold degenerate sites for estimating the neutral model

######################################################################

for chr in {1..29} Z W 
do    
    if [ -e ${pathAlignments}/${alnprefix}_${chr}.maf ]; then
	
	if [ -e ${pathResults}/${alnprefix}_${chr}_4d-codons.ss ]&&[ -e ${pathResults}/${alnprefix}_${chr}_4d-sites.ss ]; then
	    echo "already done"
	else
	    # extract codons containing 4-fold degenerate sites 
	    
    	    msa_view ${pathAlignments}/${alnprefix}_${chr}.maf --4d --features ${pathAnnot}/${annotprefix}_${chr}.gtf > ${pathResults}/${alnprefix}_${chr}_4d-codons.ss
	    
            # then extract only 4-fold degenerate sites
	    
	    msa_view ${pathResults}/${alnprefix}_${chr}_4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > ${pathResults}/${alnprefix}_${chr}_4d-sites.ss 
	    
	fi
    fi
done

######################################################################

