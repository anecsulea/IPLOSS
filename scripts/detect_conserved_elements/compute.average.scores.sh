#!/bin/bash

export dataset=$1
export coverage=$2
export cluster=$3

######################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/IPLOSS
fi

export pathResults=${path}/results/conserved_elements/${dataset}/phastCons_coverage${coverage}
export pathScripts=${path}/scripts/detect_conserved_elements

######################################################################

for chr in {1..29} Z W
do
    if [ -e ${pathResults}/most_conserved_${chr}.txt ]; then
	echo "already done"
    else
	if [ -e ${pathResults}/most_conserved_${chr}.bed ]; then
	    perl ${pathScripts}/compute.average.scores.pl --chr=${chr} --pathScores=${pathResults}/phastcons_scores_${chr}.wig --pathCoords=${pathResults}/most_conserved_${chr}.bed --pathOutput=${pathResults}/most_conserved_${chr}.txt
	else
	    echo "cannot find conserved elements for "${chr}
	fi
    fi
done

######################################################################
