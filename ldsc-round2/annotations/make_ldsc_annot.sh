#!/bin/bash

DIR=$1
PREFIX=$2

# Get annot file per annotation
header=""
for bed in $(ls ${DIR}/beds/*)
do
    header="${header},$(basename $bed .bed)"
    for chrom in $(seq 1 22)
    do
	cmd="python /storage/resources/source/ldsc/make_annot.py \
	    --bed-file $bed \
	    --bimfile /storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
	    --annot-file ${DIR}/annot_raw/$(basename $bed .bed).${chrom}.annot.gz"
	echo $cmd
    done
done | xargs -n1 -I% -P10 sh -c "%"

header=$(echo $header | sed 's/^,//')

# Merge per chromosome
for chrom in $(seq 1 22)
do
    outf=${DIR}/${PREFIX}.${chrom}.annot
    echo $header | sed 's/,/\t/g' > ${outf}
    command=paste
    for i in ${DIR}/annot_raw/*.${chrom}.annot.gz; do
	command="$command <(gzip -cd $i | cut -f5-6)"
    done
    eval $command | grep -v ANNOT >> ${outf}
    gzip -f ${outf}
done
