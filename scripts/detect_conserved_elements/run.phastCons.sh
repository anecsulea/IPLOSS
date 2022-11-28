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

for chrmaf in {1..8}
do
    export chrtype="macro_chromosomes"
    
    if [ -e ${pathResults}/most_conserved_${chrmaf}.bed ]; then
	echo "### ${chrmaf} already done ! ###"
    else
	phastCons --target-coverage 0.25 --expected-length 50 --rho 0.2 --estimate-rho ${pathResults}/rho_${chrmaf} --most-conserved ${pathResults}/most_conserved_${chrmaf}.bed --score ${pathAln}/${alnprefix}_${chrmaf}.maf ${pathMod}/phyloFit_nonconserved_4d-sites_avg_${chrtype}.mod > ${pathResults}/phastcons_scores_${chrmaf}.wig 
    fi
done

######################################################################

for chrmaf in {9..16} {18..29}
do
    export chrtype="micro_chromosomes"
    
    if [ -e ${pathResults}/most_conserved_${chrmaf}.bed ]; then
	echo "### ${chrmaf} already done ! ###"
    else
	phastCons --target-coverage 0.25 --expected-length 50 --rho 0.2 --estimate-rho ${pathResults}/rho_${chrmaf} --most-conserved ${pathResults}/most_conserved_${chrmaf}.bed --score ${pathAln}/${alnprefix}_${chrmaf}.maf ${pathMod}/phyloFit_nonconserved_4d-sites_avg_${chrtype}.mod > ${pathResults}/phastcons_scores_${chrmaf}.wig 
    fi
done

######################################################################

for chrmaf in Z
do
    export chrtype="sex_chromosomes"
    
    if [ -e ${pathResults}/most_conserved_${chrmaf}.bed ]; then
	echo "### ${chrmaf} already done ! ###"
    else
	phastCons --target-coverage 0.25 --expected-length 50 --rho 0.2 --estimate-rho ${pathResults}/rho_${chrmaf} --most-conserved ${pathResults}/most_conserved_${chrmaf}.bed --score ${pathAln}/${alnprefix}_${chrmaf}.maf ${pathMod}/phyloFit_nonconserved_4d-sites_avg_${chrtype}.mod > ${pathResults}/phastcons_scores_${chrmaf}.wig 
    fi
done

######################################################################
