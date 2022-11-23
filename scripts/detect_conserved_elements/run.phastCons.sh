#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export dataset=$1
export cluster=$2

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathAln=${path}/results/whole_genome_alignments/${dataset}/ 
export pathMod=${path}/results/conserved_elements/${dataset}/mod
export pathResults=${path}/results/conserved_elements/${dataset}/phastCons

if [ ! -d ${pathResults} ]; then
    mkdir -p ${pathResults}
fi

export alnprefix="aln"

######################################################################

for chrtype in macro_chromosomes micro_chromosomes sex_chromosomes all_chromosomes
do
    export chrmaf=${chrtype}

    if [ ${chrtype} = "all_chromosomes" ]; then
	export chrmaf="scaffolds"
    fi
    
    if [ ! -e ${pathMod}/phyloFit_nonconserved_4d-sites_avg_${chrtype}.mod ]; then
	echo  ${pathMod}/phyloFit_nonconserved_4d-sites_avg_${chrtype}.mod "is missing"
	exit
    fi
    
    if [ -e ${pathResults}/most_conserved_${chrmaf}.bed ]; then
	echo "### ${chrmaf} already done ! ###"
    else
	phastCons --target-coverage 0.25 --expected-length 50 --rho 0.2 --estimate-rho ${pathResults}/rho_${chrmaf} --most-conserved ${pathResults}/most_conserved_${chrmaf}.bed --score ${pathAln}/${alnprefix}_${chrmaf}.maf ${pathMod}/phyloFit_nonconserved_4d-sites_avg_${chrtype}.mod > ${pathResults}/phastcons_scores_${chrmaf}.wig 
    fi
done

######################################################################
