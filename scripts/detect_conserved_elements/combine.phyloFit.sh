#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export dataset=$1
export cluster=$2
export nthreads=$3

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathResults=${path}/results/conserved_elements/${dataset}/phyloFit/

######################################################################

## macro-chromosomes

export pathsNonConserved=""

for chr in {1..8}
do
    export pathsNonConserved="${pathResults}/phyloFit_nonconserved_${chr}_4d-sites.mod",${pathsNonConserved} 
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_macro_chromosomes.mod

######################################################################

## micro-chromosomes

export pathsNonConserved=""

for chr in {9..16} {18..29} 
do
    export pathsNonConserved="${pathResults}/phyloFit_nonconserved_${chr}_4d-sites.mod",${pathsNonConserved} 
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_micro_chromosomes.mod

######################################################################

## sex-chromosomes (simple cp for phast pipeline homogenisation)

cp ${pathResults}/phyloFit_nonconserved_chrZ_4d-sites.mod ${pathResults}/phyloFit_nonconserved_4d-sites_avg_sex_chromosomes.mod

######################################################################

## micro/macro/sex (for scaffold)

export pathsNonConserved=""

for chr in {1..16} {18..29} Z
do
    export pathsNonConserved="${pathResults}/phyloFit_nonconserved_${chr}_4d-sites.mod",${pathsNonConserved} 
done

phyloBoot --read-mods ${pathsNonConserved} --output-average ${pathResults}/phyloFit_nonconserved_4d-sites_avg_all_chromosomes.mod

######################################################################


