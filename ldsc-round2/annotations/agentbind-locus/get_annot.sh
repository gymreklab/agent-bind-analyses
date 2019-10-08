#!/bin/bash

factor=$1
top1=$2
top2=$3
top3=$4
top4=$5

OUTDIR=/storage/mgymrek/agent-bind/ldsc/annotations/agentbind-locus/beds

rm ${OUTDIR}/${factor}.*

num=1
for score in $top1 $top2 $top3 $top4
do
    cat /storage/pandaman/project/AgentBind-LD/results/c/${factor}/lable-logit.txt | \
	awk -v"score=$score" '($NF>=score) {print $1 "\t" $2-1 "\t" $2+1}' | awk '($2>0)' | awk '{print "chr"$0}' > \
	${OUTDIR}/${factor}.top${num}.bed
    num=$((num+1))
done

