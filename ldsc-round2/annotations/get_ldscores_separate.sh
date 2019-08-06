#!/bin/bash

DIR=$1
for bed in $(ls ${DIR}/beds/*)
do
    for chrom in $(seq 1 22)
    do
	cmd="python /storage/resources/source/ldsc/ldsc.py \
	    --l2 \
	    --bfile /storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
	    --ld-wind-cm 1 \
	    --annot ${DIR}/annot_raw/$(basename $bed .bed).${chrom}.annot.gz \
	    --out ${DIR}/annot_raw/$(basename $bed .bed).${chrom} \
	    --thin-annot \
	    --print-snps /storage/mgymrek/agent-bind/ldsc/resources/hapmap3_snps_v2/hm.${chrom}.snp"
	echo $cmd
    done
done | xargs -I% -n1 -P10 sh -c "%"
