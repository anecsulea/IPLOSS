#!/bin/bash

export sp1=$1
export sp2=$2
export annot=$3
export cluster=$4
export modulo=$5
export run=$6

####################################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/IPLOSS
fi

echo ${cluster} ${path}

####################################################################################

export pathProjections=${path}/results/liftOver_gene_families
export pathStringTie=${path}/results/stringtie_assembly
export pathEnsembl=${path}/data/ensembl_annotations
export pathResults=${path}/results/liftOver_gene_families
export pathScripts=${path}/scripts/gene_families

export release=103

####################################################################################

if [ ${annot} = "Ensembl" ]; then
    export prefix=ExonBlocks_FilteredTranscripts_Ensembl${release}
    export pathExons=${pathEnsembl}
fi

if [ ${annot} = "EnsemblStringTie" ]; then
    export prefix=ExonBlocks_combined_annotations_StringTie_Ensembl
    export pathExons=${pathStringTie}
fi

####################################################################################

if [ ${sp1} = ${sp2} ]; then
    echo "cannot project from "${sp1}" to "${sp2}
    exit
fi

####################################################################################

if [ -e ${pathResults}/pecan_alignments_projection_clusters/${sp1}_${sp2}_${annot} ]; then
    echo "dir output already there"
else
    mkdir -p ${pathResults}/pecan_alignments_projection_clusters/${sp1}_${sp2}_${annot}
fi

####################################################################################

if [ ${cluster} = "in2p3" ]; then
    for i in {0..49}
    do
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_pecan
	echo "#SBATCH --job-name=pecan_${i}" >>  ${pathScripts}/bsub_script_pecan
	echo "#SBATCH --output=${pathScripts}/std_output_pecan_${i}.txt" >>  ${pathScripts}/bsub_script_pecan
	echo "#SBATCH --error=${pathScripts}/std_error_pecan_${i}.txt" >> ${pathScripts}/bsub_script_pecan
	echo "#SBATCH --ntasks=1" >> ${pathScripts}/bsub_script_pecan
	echo "#SBATCH --cpus-per-task=1" >> ${pathScripts}/bsub_script_pecan
	echo "#SBATCH --time=7-00:00:00" >> ${pathScripts}/bsub_script_pecan

	echo "perl ${pathScripts}/run.pecan.projection.clusters.pl --species1=${sp1} --species2=${sp2} --pathSequences1=${pathExons}/${sp1}/${prefix}_cDNASequences.fa --pathSequences2=${pathExons}/${sp2}/${prefix}_cDNASequences.fa --pathClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_${prefix}.txt --run=${i} --modulo=50 --dirOutput=${pathResults}/pecan_alignments_projection_clusters/${sp1}_${sp2}_${annot}" >> ${pathScripts}/bsub_script_pecan

	sbatch ${pathScripts}/bsub_script_pecan

    done
fi

####################################################################################

if [ ${cluster} = "cloud" ]; then
    perl ${pathScripts}/run.pecan.projection.clusters.pl --species1=${sp1} --species2=${sp2} --pathSequences1=${pathExons}/${sp1}/${prefix}_cDNASequences.fa --pathSequences2=${pathExons}/${sp2}/${prefix}_cDNASequences.fa --pathClusters=${pathResults}/ProjectionClusters_${sp1}_${sp2}_${prefix}.txt --run=${run} --modulo=${modulo} --dirOutput=${pathResults}/pecan_alignments_projection_clusters/${sp1}_${sp2}_${annot}
fi

####################################################################################
