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
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
    export pathTools=/ifb/data/mydatalocal/Tools
fi

export pathResults=${path}/results/CDS_annotation/${sp}/GeMoMa

#########################################################################

if [ -e ${pathResults}/combined ]; then
    echo "results dir already there"
else
    mkdir -p ${pathResults}/combined
fi

#########################################################################

echo "#!/bin/bash " > script_combine_GeMoMa
echo -n "java -jar ${pathTools}/GeMoMa/GeMoMa-1.9.jar CLI GAF " >> script_combine_GeMoMa

#########################################################################

for ref in `ls ${pathResults} | grep -v combined`
do
    if [ -e ${pathResults}/${ref}/final_annotation.gff ]; then
	echo -n "p=${ref} g=${pathResults}/${ref}/final_annotation.gff ">> script_combine_GeMoMa
    fi
done

#########################################################################

echo  "outdir=${pathResults}/combined " >> script_combine_GeMoMa

#########################################################################

chmod a+x script_combine_GeMoMa

./script_combine_GeMoMa

#########################################################################
