#!/bin/bash

export sp=$1
export cluster=$2

#########################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/${USER}/IPLOSS
    export pathTools=/beegfs/home/${USER}/IPLOSS
fi

if [ ${cluster} = "in2p3" ]; then
    export path=/sps/biometr/necsulea/IPLOSS
    export pathTools=/sps/biometr/necsulea/Tools
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
    export pathTools=/ifb/data/mydatalocal/Tools
fi

export release=103
export pathGenomeSequence=${path}/data/genome_sequences/${sp}/genome_ensembl${release}.fa
export pathGeMoMa=${path}/results/CDS_annotation/${sp}/GeMoMa/combined
export pathStringTie=${path}/results/stringtie_assembly/${sp}
export pathProteinSequences=${path}/data/protein_sequences/Ensembl${release}
export pathScripts=${path}/scripts/CDS_annotation

#########################################################################

if [ -e ${pathGeMoMa}/filtered_predictions.gtf ]; then
    echo "GTF file already there"
else
    gffread ${pathGeMoMa}/filtered_predictions.gff -T -o  ${pathGeMoMa}/filtered_predictions.gtf 
fi

#########################################################################

if [ -e ${pathGeMoMa}/filtered_predictions.faa ]; then
    echo "Fasta file for proteins already there"
else
    gffread ${pathGeMoMa}/filtered_predictions.gff -g ${pathGenomeSequence} -y ${pathGeMoMa}/filtered_predictions.faa
fi

#########################################################################

if [ ${sp} == "Chicken" ]; then
    export pathKnownProteins=${pathProteinSequences}/Gallus_gallus.GRCg6a.pep.all.fa
fi

if [ ${sp} == "Duck" ]; then
    export pathKnownProteins=${pathProteinSequences}/Anas_platyrhynchos_platyrhynchos.CAU_duck1.0.pep.all.fa 
fi

#########################################################################

perl ${pathScripts}/extract.compatible.CDS.pl --pathNewCDSAnnotation=${pathGeMoMa}/filtered_predictions.gtf --pathTranscriptAnnotation=${pathStringTie}/combined_annotations_StringTie_Ensembl.gtf --minFractionOverlap=0.5 --pathKnownProteins=${pathKnownProteins} --pathNewProteins=${pathGeMoMa}/filtered_predictions.faa --pathOutputProteins=${pathGeMoMa}/combined_protein_sequences.faa --pathOutputOverlap=${pathGeMoMa}/compatible_CDS.txt

#########################################################################
