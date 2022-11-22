#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export refsp=$1
export dataset=$2
export cluster=$3

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathAln=${path}/results/whole_genome_alignments/${dataset}
export pathSites=${path}/results/conserved_elements/${dataset}/4d/
export pathResults=${path}/results/conserved_elements/${dataset}/mod/

######################################################################

for chr in {1..29} Z W 
do
    export alnprefix="aln"
    
    if [ -e ${pathResults}/phyloFit_nonconserved_4d-sites_${chr}.mod ]; then
	echo "already done"
    else
	phyloFit --tree ${pathAln}/tree.txt --msa-format SS --out-root ${pathResults}/phyloFit_nonconserved_4d-sites_${chr} ${pathSites}/${alnprefix}_${chr}_4d-sites.ss
    fi
done

######################################################################
