#!/bin/bash

source params.sh
CHROMS=$(seq 1 22)

top1=0.43 ;top5=0.1256; topall=0 # TODO get automatically

#for chrom in $CHROMS
#do
#    echo "CHR,BP,SNP,CM,IMPACT-top1,IMPACT-top5,IMPACT-topall" | sed 's/,/\t/g' > ${OUTDIR}/annot/IMPACT.${chrom}.annot
#    zcat /storage/cynthiawu/soft/ldsc/RA_Analysis/Annot/IMPACT/Annot_HapMapOnly/IMPACT_HapMap_chr.${chrom}.annot.gz | grep -v CHR | awk -v"top1=$top1" -v"top5=$top5" -v"topall=$topall" \
#	'{print $1 "\t" $2 "\t" $3 "\t0\t" ($5>top1?"1":"0") "\t" ($5>top5?"1":"0") "\t" ($5>topall?1:0)}' >># ${OUTDIR}/annot/IMPACT.${chrom}.annot
#    gzip ${OUTDIR}/annot/IMPACT.${chrom}.annot
#done

for chrom in $CHROMS
do
    echo "python /storage/cynthiawu/soft/ldsc/ldsc.py \
	--l2 \
	--bfile ${BFILE}.${chrom} \
	--ld-wind-cm 1 \
	--annot ${OUTDIR}/annot/IMPACT.${chrom}.annot.gz \
	--out ${OUTDIR}/annot/IMPACT.${chrom} --print-snps ${HMAPSNPS}"
done | xargs -n1 -I% -P22 sh -c "%"
