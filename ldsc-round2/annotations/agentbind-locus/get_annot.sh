#!/bin/bash

factor=$1
top5=$2
top1=$3
OUTDIR=/storage/mgymrek/agent-bind/ldsc/annotations/agentbind-locus/beds

cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | \
    awk -v"score=$top5" '($NF>=score) {print $1 "\t" $2-500 "\t" $2+500}' | awk '($2>0)' > \
    ${OUTDIR}/${factor}.top5.bed

cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | \
    awk -v"score=$top1" '($NF>=score) {print $1 "\t" $2-500 "\t" $2+500}' | awk '($2>0)' > \
    ${OUTDIR}/${factor}.top1.bed
