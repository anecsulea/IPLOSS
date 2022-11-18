#!/bin/bash

export target=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
    export pathTools=/beegfs/home/${USER}/Tools
    export version=1.9
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/IPLOSS
    export pathTools=/sps/biometr/necsulea/Tools
    export version=1.9
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools
    export version=1.9
fi

export release=103

export pathAllGenomes=${path}/data/genome_sequences
export pathGenomeSequence=${path}/data/genome_sequences/${target}/genome_ensembl${release}.fa
export pathRNASeq=${path}/results/RNASeq_alignments/${target}/all_samples
export pathResults=${path}/results/CDS_annotation/${target}/GeMoMa
export pathScripts=${path}/scripts/CDS_annotation

## using mmseqs2, installed using apt on Ubuntu 20.04
## version 9-d36de+ds-4 amd64

#########################################################################

if [ -e ${pathResults} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}
fi

#########################################################################

echo "#!/bin/bash" > ${pathScripts}/bsub_script_gemoma

#############################################

if [ ${cluster} = "pbil" ]; then
    echo "#SBATCH --job-name=gemoma_rna_seq_${target}" >>  ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${target}.txt" >>  ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${target}.txt" >> ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_gemoma
    
    echo "singularity exec -B ${path} -B ${pathTools} ${pathTools}/basic_ubuntu.simg java -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI ERE s=FR_FIRST_STRAND m=${pathRNASeq}/accepted_hits.bam outdir=${pathResults}" >> ${pathScripts}/bsub_script_gemoma
    
    sbatch ${pathScripts}/bsub_script_gemoma
fi

#############################################

if [ ${cluster} = "in2p3" ]; then
    echo "#SBATCH --job-name=gemoma_${target}" >>  ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${target}.txt" >>  ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${target}.txt" >> ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_gemoma
    echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_gemoma
    
    echo "java -Xms2G -Xmx64G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI ERE s=FR_FIRST_STRAND m=${pathRNASeq}/accepted_hits.bam outdir=${pathResults}" >> ${pathScripts}/bsub_script_gemoma
    
    sbatch ${pathScripts}/bsub_script_gemoma
fi

#############################################

if [ ${cluster} = "cloud" ]; then
    ## mmseqs available in PATH
    echo "java  -Xms2G -Xmx64G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI ERE s=FR_FIRST_STRAND m=${pathRNASeq}/accepted_hits.bam outdir=${pathResults}" >> ${pathScripts}/bsub_script_gemoma
    
    chmod a+x ${pathScripts}/bsub_script_gemoma
    ${pathScripts}/bsub_script_gemoma
fi

#########################################################################
