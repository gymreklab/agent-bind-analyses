#!/bin/bash

factor=$1
top5=$2
top1=$3

OUTDIR=/storage/mgymrek/agent-bind/ldsc/annotations/agentbind-locus/beds

cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | \
    awk -v"score=$top5" '($NF>=score) {print $1 "\t" $2-1 "\t" $2+1}' | awk '($2>0)' | awk '{print "chr"$0}' > \
    ${OUTDIR}/${factor}.top5.bed

cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | \
    awk -v"score=$top1" '($NF>=score) {print $1 "\t" $2-1 "\t" $2+1}' | awk '($2>0)' | awk '{print "chr"$0}' > \
    ${OUTDIR}/${factor}.top1.bed

cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | \
    awk '($NF>0.001) {print $1 "\t" $2-1 "\t" $2+1}' | awk '($2>0)' | awk '{print "chr"$0}' > \
    ${OUTDIR}/${factor}.all.bed

cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | shuf | head -n 50000 | \
    awk '{print $1 "\t" $2-1 "\t" $2+1}' | awk '($2>0)' | awk '{print "chr"$0}' > \
    ${OUTDIR}/${factor}.random.bed
