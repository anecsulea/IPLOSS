#!/bin/bash

## original script by Alexandre Laverr� & Anamaria Necsulea

export cluster=$1
export nthreads=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathHAL=${path}/results/whole_genome_alignments/avian_366
export pathResults=${path}/results/whole_genome_alignments/species_with_IP

#########################################################################

export targetGenomes="Chauna_torquata,Anser_cygnoid,Cairina_moschata,Asarcornis_scutulata,Anas_zonorhyncha,Anas_platyrhynchos,Anas_platyrhynchos_platyrhynchos,Anseranas_semipalmata,Pauxi_pauxi,Penelope_pileata,Casuarius_casuarius,Dromaius_novaehollandiae,Eudromia_elegans,Nothocercus_nigrocapillus,Nothocercus_julius,Nothoprocta_perdicaria,Nothoprocta_ornata,Nothoprocta_pentlandii,Tinamus_guttatus,Rhea_americana,Rhea_pennata,Struthio_camelus,"

#########################################################################

docker run -v ${path}:/ifb/data/mydatalocal/IPLOSS --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 hal2mafMP.py ${pathHAL}/366-avian.hal ${pathResults}/aln.maf  --targetGenomes ${targetGenomes} --refGenome Anas_platyrhynchos_platyrhynchos --numProc ${nthreads} --noDupes --splitBySequence --smallSize 100000

#########################################################################
