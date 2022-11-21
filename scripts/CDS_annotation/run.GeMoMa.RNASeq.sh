#!/bin/bash

export target=$1
export ref=$2
export source=$3
export cluster=$4
export threads=$5
export part=$6

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

#########################################################################

export release=103

export pathAllGenomes=${path}/data/genome_sequences
export pathGenomeSequence=${path}/data/genome_sequences/${target}/


if [ ${source} = "Reptiles_Ensembl103" ]; then
    export pathSourceAnnotations=${path}/data/ensembl_annotations/${source}
    export pathSourceGenomes=${path}/data/genome_sequences/${source}
fi

if [ ${source} = "NCBI" ]; then
    export pathSourceAnnotations=${path}/data/NCBI_annotations
    export pathSourceGenomes=${path}/data/genome_sequences/${source}
fi


if [ ${source} = "Reptiles_Ensembl103_parts" ]; then
    export pathSourceAnnotations=${path}/data/ensembl_annotations/Reptiles_Ensembl103/parts
    export pathSourceGenomes=${path}/data/genome_sequences/Reptiles_Ensembl103
fi

if [ ${source} = "NCBI_parts" ]; then
    export pathSourceAnnotations=${path}/data/NCBI_annotations/parts
    export pathSourceGenomes=${path}/data/genome_sequences/NCBI
fi


export pathResults=${path}/results/CDS_annotation/${target}/GeMoMa
export pathScripts=${path}/scripts/CDS_annotation

## in2p3 MMseqs2 Version: 9cc89aa594131293b8bc2e7a121e2ed412f0b931

#########################################################################

if [ -e ${pathGenomeSequence}/genome_ensembl${release}.clean.fa ]; then
    echo "clean target genome sequence already there"
else
    perl ${pathScripts}/cleanup.fasta.names.pl --pathInput=${pathGenomeSequence}/genome_ensembl${release}.fa --pathOutput=${pathGenomeSequence}/genome_ensembl${release}.clean.fa 
fi

export pathTargetGenome=${pathGenomeSequence}/genome_ensembl${release}.clean.fa

echo "target genome" ${pathTargetGenome}

#########################################################################

export refGenome=`ls ${pathSourceGenomes} | grep ${ref} | grep .fa.gz`
export prefix=`basename ${refGenome} .fa.gz`

if [ -e ${pathSourceGenomes}/${prefix}.clean.fa ];then
    echo "clean reference genome sequence already there"
else
     perl ${pathScripts}/cleanup.fasta.names.pl --pathInput=${pathSourceGenomes}/${prefix}.fa.gz --pathOutput=${pathSourceGenomes}/${prefix}.clean.fa 
fi

export pathRefGenome=${pathSourceGenomes}/${prefix}.clean.fa 

echo "reference genome" ${pathRefGenome}

#########################################################################

if [ ${source} = "Reptiles_Ensembl103" ]||[ ${source} = "NCBI" ]; then
    export annotfile=`ls ${pathSourceAnnotations} | grep ${ref}'\.' | grep gff`
    export pathRefAnnot=${pathSourceAnnotations}/${annotfile}
    export outdir=${ref}
fi


if [ ${source} = "Reptiles_Ensembl103_parts" ]||[ ${source} = "NCBI_parts" ]; then
    export annotfile=`ls ${pathSourceAnnotations} | grep ${ref}'\.' | grep part${part}.gff`
    export pathRefAnnot=${pathSourceAnnotations}/${annotfile}
    export outdir=${ref}_part${part}
fi

echo "reference annotation" ${pathRefAnnot}

if [ -e ${pathRefAnnot} ]; then
    echo "ok, path annotation exists"
else
    echo "cannot find reference annotation"
    exit
fi

#########################################################################

if [ -e ${pathResults}/${outdir} ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}/${outdir}
fi

#########################################################################

if [ -e ${pathResults}/${outdir}/final_annotation.gff ]; then
    echo "already done"
else
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_gemoma

    #############################################

    if [ ${cluster} = "pbil" ]; then
	echo "#SBATCH --job-name=gemoma_${ref}" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${ref}_${target}_part${part}.txt" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${ref}_${target}_part${part}.txt" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --partition=normal" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --mem=12G" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --time=24:00:00" >> ${pathScripts}/bsub_script_gemoma

	echo "singularity exec -B ${path} -B ${pathTools} ${pathTools}/basic_ubuntu.simg java -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults}/${outdir} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathTargetGenome} i=${ref} a=${pathRefAnnot}  g=${pathRefGenome} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 m=${pathTools}/mmseqs/bin/ r=EXTRACTED introns=${pathResults}/introns.gff coverage=STRANDED coverage_forward=${pathResults}/coverage_forward.bedgraph coverage_reverse=${pathResults}/coverage_reverse.bedgraph " >> ${pathScripts}/bsub_script_gemoma

	sbatch ${pathScripts}/bsub_script_gemoma
    fi

    #############################################

    if [ ${cluster} = "in2p3" ]; then
	echo "#SBATCH --job-name=gemoma_${ref}" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --output=${pathScripts}/std_output_GEMOMA_${ref}_${target}_part${part}.txt" >>  ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --error=${pathScripts}/std_error_GEMOMA_${ref}_${target}_part${part}.txt" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --cpus-per-task=${threads}" >> ${pathScripts}/bsub_script_gemoma
	echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_gemoma

	echo "java -Xms2G -Xmx64G  -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults}/${outdir} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathTargetGenome} i=${ref} a=${pathRefAnnot}  g=${pathRefGenome} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10 m=${pathTools}/mmseqs/bin/ r=EXTRACTED introns=${pathResults}/introns.gff coverage=STRANDED coverage_forward=${pathResults}/coverage_forward.bedgraph coverage_reverse=${pathResults}/coverage_reverse.bedgraph " >> ${pathScripts}/bsub_script_gemoma

	sbatch ${pathScripts}/bsub_script_gemoma
    fi


    #############################################

    if [ ${cluster} = "cloud" ]; then
	echo "java  -Xms2G -Xmx64G -jar ${pathTools}/GeMoMa/GeMoMa-${version}.jar CLI GeMoMaPipeline threads=${threads} outdir=${pathResults}/${outdir} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=${pathTargetGenome} i=${ref} a=${pathRefAnnot}  g=${pathRefGenome} GeMoMa.m=500000 Extractor.f=false GeMoMa.i=10  r=EXTRACTED introns=${pathResults}/introns.gff coverage=STRANDED coverage_forward=${pathResults}/coverage_forward.bedgraph coverage_reverse=${pathResults}/coverage_reverse.bedgraph " >> ${pathScripts}/bsub_script_gemoma

	chmod a+x ${pathScripts}/bsub_script_gemoma
	${pathScripts}/bsub_script_gemoma
    fi

fi

#########################################################################
