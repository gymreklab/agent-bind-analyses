#!/bin/bash

set -e

DIR=/storage/mgymrek/agent-bind/ldsc/annotations/agentbind-locus-quant/
# Get SNP lists
#for chrom in $(seq 1 22)
#do
#    cat /storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim  | awk '{print $1 "\t" $4 "\t" $4+1}' > ${DIR}/snps/${chrom}.bed
#done

# Make annotation files

for factor in Gata3 Foxp3 Stat3
do
    for chrom in $(seq 1 22)
    do
	echo $chrom
	outfile=${DIR}/annot_raw/${factor}.${chrom}.annot
#	./get_annot.py /storage/pandaman/project/AgentBind-LD/results-full/c/${factor}/lable-logit.txt ${chrom} /storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim > ${outfile}
#	./get_annot.py /storage/pandaman/project/AgentBind-LD/results-full/b/${factor}/lable-logit.txt ${chrom} /storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim > ${outfile}
    done
done

# Merge to a single file per chrom
PREFIX=AgentBindLocusQuant
header="Gata3,Foxp3,Stat3"
for chrom in $(seq 1 22)
do
    outf=${DIR}/${PREFIX}.${chrom}.annot
    echo $header | sed 's/,/\t/g' > ${outf}
    f1=${DIR}/annot_raw/Gata3.${chrom}.annot
    f2=${DIR}/annot_raw/Foxp3.${chrom}.annot
    f3=${DIR}/annot_raw/Stat3.${chrom}.annot
    paste $f1 $f2 $f3 | grep -v ANNOT >> ${outf}
    gzip -f ${outf}
done

# Compute LD scores
../get_ldscores.sh ${DIR} AgentBindLocusQuant
