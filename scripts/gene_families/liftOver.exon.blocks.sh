#!/bin/bash

export ref=$1
export tg=$2
export annot=$3
export cluster=$4

######################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

######################################################################

export pathEnsembl=${path}/data/ensembl_annotations/${ref}
export pathStringTie=${path}/results/stringtie_assembly/${ref}
export pathGenomeAlignments=${path}/results/whole_genome_alignments/pairwise_Gallus_gallus_vs_Anas_platyrhynchos_platyrhynchos
export pathResults=${path}/results/liftOver_gene_families
export pathScripts=${path}/scripts/gene_families

export release=103

######################################################################

if [ ${annot} = "EnsemblStringTie" ]; then
    export pathExons=${pathStringTie}
    export prefixExons=ExonBlocks_combined_annotations_StringTie_Ensembl
fi

######################################################################

if [ ${ref} = "Chicken" ]; then
    export outref="Gallus_gallus"
fi

if [ ${ref} = "Duck" ]; then
    export outref="Anas_platyrhynchos_platyrhynchos"
fi

######################################################################

if [ ${tg} = "Chicken" ]; then
    export outtg="Gallus_gallus"
fi

if [ ${tg} = "Duck" ]; then
    export outtg="Anas_platyrhynchos_platyrhynchos"
fi

######################################################################

if [ -e ${pathResults} ]; then
    echo "path results already there"
else
    mkdir -p ${pathResults}
fi

######################################################################

for file in `ls ${pathExons} | grep ${prefixExons} | grep bed`
do
    export prefix=`basename ${file} .bed`
    echo ${file} ${prefix}

    if [ -e  ${pathResults}/From${ref}_To${tg}_${prefix}.bed ]; then
	echo "already done"
    else
	echo "#!/bin/bash" >  ${pathScripts}/bsub_script_liftOver

	echo "docker run -v ${path}:${path} --rm -t quay.io/comparative-genomics-toolkit/cactus:v2.0.3 halLiftover ${pathGenomeAlignments}/alignment.hal ${outref} ${pathExons}/${prefix}.bed ${outtg} ${pathResults}/From${ref}_To${tg}_${prefix}.bed " >> ${pathScripts}/bsub_script_liftOver

	if [ ${cluster} = "cloud" ]; then
	    chmod a+x ${pathScripts}/bsub_script_liftOver
	    ${pathScripts}/bsub_script_liftOver
	fi
    fi
done

######################################################################
