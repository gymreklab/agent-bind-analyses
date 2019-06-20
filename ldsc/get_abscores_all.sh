#!?bin/bash

source params.sh
CHROMS=$(seq 1 22)

for TF in Gata3 Foxp3 #GATA3 FOXP3 STAT3 TBET
do
    for chrom in $CHROMS
    do
	# Get annot files with top % of scores - AgentBind Importance
	echo "./get_abscores.sh $chrom $TF ${TF^^}"
    done 
done | xargs -n1 -I % -P22 sh -c "%"

for TF in Stat3 Tbet
do
    for chrom in $CHROMS
    do
	# Get annot files with top % of scores - AgentBind Importance
	echo "./get_abscores.sh $chrom $TF $TF"
    done 
done | xargs -n1 -I % -P22 sh -c "%"
