#!/bin/bash

## adapted from script by Menghan Wang

export sp=$1
export cluster=$2

####################################################################################

if [ ${cluster} = "pbil" ]||[ ${cluster} = "pbillocal" ]; then
    export path=/beegfs/data/necsulea/IPLOSS
fi

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathGenomes=${path}/data/genome_sequences/${sp}
export pathResults=${path}/results/snATACSeq_indexes/${sp}
export pathCDSannotation=${path}/scripts/CDS_annotation

export release=103

####################################################################################

## keep only assembled chromosomes

if [ -e ${pathResults}/chromosomes.txt ]; then
    echo "chromosome list already there, not doing anything"
else
    
    if [ ${sp} = "Chicken" ]; then
	export chrlist="c(\"1\""
	echo 1 >> ${pathResults}/chromosomes.txt
	
	for chr in {2..28} {30..33} MT Z W
	do
	    echo ${chr} >> ${pathResults}/chromosomes.txt
	    export chrlist=${chrlist}", \""${chr}"\""
	done

	export chrlist=${chrlist}")"
    fi
    
    if [ ${sp} = "Duck" ]; then
	export chrlist="c(\"1\""
	echo 1 >> ${pathResults}/chromosomes.txt

	for chr in {2..16} {18..29} Z
	do	  
	    echo ${chr} >> ${pathResults}/chromosomes.txt
	    export chrlist=${chrlist}", \""${chr}"\""
	done
	
	export chrlist=${chrlist}")"
    fi
fi

####################################################################################

if [ -e ${pathGenomes}/genome_ensembl${release}.clean.fa ]; then
    echo "clean fasta file already there"
else
    perl ${pathCDSannotation}/cleanup.fasta.pl --pathInput=${pathGenomes}/genome_ensembl${release}.fa --pathOutput=${pathGenomes}/genome_ensembl${release}.clean.fa
fi

####################################################################################

if [ -e ${pathResults}/genome_ensembl${release}_chromosomes.fa ]; then
    echo "filtered chromosomes already there, not doing anything"
else
    seqfilter -i ${pathGenomes}/genome_ensembl${release}.clean.fa -l ${pathResults}/chromosomes.txt -o ${pathResults}/genome_ensembl${release}_chromosomes.fa
fi

####################################################################################

if [ -e ${pathResults}/genome_ensembl${release}_chromosomes.2bit ]; then
    echo "2bit file already there, not doing anything"
else
    faToTwoBit ${pathResults}/genome_ensembl${release}_chromosomes.fa ${pathResults}/genome_ensembl${release}_chromosomes.2bit
fi

####################################################################################

if [ ${sp} = "Chicken" ]; then
    export latinname="Gallus gallus"
    export abbr="Ggallus"
    export biocview="Gallus_gallus"
fi

if [ ${sp} = "Duck" ]; then
    export latinname="Anas platyrhynchos platyrhynchos"
    export abbr="Aplatyrhynchos"
    export biocview="Anas_platyrhynchos_platyrhynchos"
fi

####################################################################################

echo "Package: BSgenome.${abbr}"  > ${pathResults}/seed_file.txt
echo "Title: Full genome of ${latinname} for scATAC" >> ${pathResults}/seed_file.txt
echo "Description: ${latinname} genome modified under 10x" >> ${pathResults}/seed_file.txt
echo "Version: 1.0.0" >> ${pathResults}/seed_file.txt
echo "organism: ${latinname}" >> ${pathResults}/seed_file.txt
echo "organism_biocview: ${biocview}" >> ${pathResults}/seed_file.txt
echo "BSgenomeObjname: ${abbr}" >> ${pathResults}/seed_file.txt
echo "Author: Menghan Wang" >> ${pathResults}/seed_file.txt
echo "common_name: ${sp}" >> ${pathResults}/seed_file.txt
echo "provider: ensembl" >> ${pathResults}/seed_file.txt

## not the actual genome for duck but we need a workaround for BSgenome
echo "provider_version: GRCg6a" >> ${pathResults}/seed_file.txt

echo "release_date: November.2022" >> ${pathResults}/seed_file.txt
echo "seqs_srcdir: ${pathResults}" >> ${pathResults}/seed_file.txt
echo "seqfile_name: genome_ensembl${release}_chromosomes.2bit" >> ${pathResults}/seed_file.txt
echo "seqnames: ${chrlist}" >> ${pathResults}/seed_file.txt

####################################################################################
