#!/bin/bash

export sp=$1
export assembly=$2
export cluster=$3
export nthreads=$4

########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

export pathIndexes=${path}/results/genome_indexes/${sp}
export pathScripts=${path}/scripts/genome_indexes

## hisat2-2.2.1

########################################################################

export prefix=genome_${assembly}
export pathGenomeSequence=${pathIndexes}/${prefix}.fa

########################################################################

echo "#!/bin/bash " > ${pathScripts}/bsub_script_build_${sp}

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=build_${sp}" >>  ${pathScripts}/bsub_script_build_${sp}
    echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_build_${sp}
    echo "#SBATCH --output=${pathScripts}/std_out_build_${sp}" >>  ${pathScripts}/bsub_script_build_${sp}
    echo "#SBATCH --error=${pathScripts}/std_err_build_${sp}" >>  ${pathScripts}/bsub_script_build_${sp}
    echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script_build_${sp}
    echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_build_${sp}
    echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_build_${sp}
fi

echo "hisat2-build --seed 19 -p ${nthreads} ${pathGenomeSequence} ${pathIndexes}/${prefix}" >> ${pathScripts}/bsub_script_build_${sp}

########################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "in2p3" ]; then
    sbatch ${pathScripts}/bsub_script_build_${sp}
fi

if [ ${cluster} = "cloud" ]; then
    chmod a+x ${pathScripts}/bsub_script_build_${sp}
    ${pathScripts}/bsub_script_build_${sp}
fi

########################################################################
