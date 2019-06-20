#!/bin/bash

set -e

chrom=$1
TF=$2
TFupper=$3 #${TF^^}

source params.sh

# Get annot files with top % of scores - AgentBind Importance
toptop=0.99; top1=0.85; top5=0.3 topall=0 # TODO get automatically for each one
echo "CHR,BP,SNP,CM,AB-top99,AB-top85,AB-top30,AB-topall" | sed 's/,/\t/g' > ${OUTDIR}/annot/${TF}.AB.${chrom}.annot
zcat /storage/cynthiawu/soft/ldsc/RA_Analysis/Annot/ImportanceScores/${TFupper}/${TFupper}chr.${chrom}.annot.gz | \
    grep -v CHR |  awk -v"toptop=$toptop" -v"top1=$top1" -v"top5=$top5" -v"topall=$topall" \
    '{print $1 "\t" $2 "\t" $3 "\t0\t" ($5>toptop?"1":"0") "\t" ($5>top1?"1":"0") "\t" ($5>top5?"1":"0") "\t" ($5>topall?1:0)}' >> ${OUTDIR}/annot/${TF}.AB.${chrom}.annot
gzip -f ${OUTDIR}/annot/${TF}.AB.${chrom}.annot
python /storage/cynthiawu/soft/ldsc/ldsc.py \
    --l2 \
    --bfile ${BFILE}.${chrom} \
    --ld-wind-cm 1 \
    --annot ${OUTDIR}/annot/${TF}.AB.${chrom}.annot.gz \
    --out ${OUTDIR}/annot/${TF}.AB.${chrom} --print-snps ${HMAPSNPS}

