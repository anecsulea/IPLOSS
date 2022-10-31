#!/bin/bash

export type=$1  ## tree inference method
export cluster=$2
export threads=$3

##########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
fi

export pathResults=${path}/results/gene_families/OrthoFinder
export pathScripts=${path}/scripts/coding_gene_families

##########################################################################

## OrthoFinder 2.5.4
## IQ-TREE multicore version 2.0.7 for Linux 64-bit 
## diamond version 2.0.5
## MUSCLE v3.8.1551

##########################################################################

if [ ${cluster} = "cloud" ]; then
    ulimit -n 100000
fi

##########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_orthofinder

##########################################################################

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=orthofinder_${type}" >>  ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --output=${pathScripts}/std_output_orthofinder.txt" >>  ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --error=${pathScripts}/std_error_orthofinder.txt" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --mem=31G" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_orthofinder
    echo "#SBATCH --time=168:00:00" >> ${pathScripts}/bsub_script_orthofinder
fi
    
##########################################################################

echo "orthofinder -f ${pathResults} -o ${pathResults}/${type} -t ${threads} -a ${threads} -I 2 -S diamond -A muscle -M msa -y -T ${type}" >>  ${pathScripts}/bsub_script_orthofinder

##########################################################################

if [ ${cluster} = "pbil" ]; then
    sbatch ${pathScripts}/bsub_script_orthofinder
fi

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_orthofinder
    ${pathScripts}/bsub_script_orthofinder
fi

##########################################################################
