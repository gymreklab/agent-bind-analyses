#!/bin/bash

factor=$1

echo $factor
mkdir -p /storage/mgymrek/agent-bind/singletons/${factor}

# Get overall distribution of raw scores
cat /storage/pandaman/project/AgentBind-GM12878-DanQ-unfixed-rnn-trans/storage/AgentBind-GM12878-DanQ/tmp/${factor}+GM12878/seqs_one_hot_c/vis-weights-total/weight.txt | awk '(NR%2==0)' | sed 's/;/\n/g' | sort -k1,1g -T /storage/mgymrek/agent-bind/singletons/${factor} > /storage/mgymrek/agent-bind/singletons/${factor}/raw_vals.txt
exit 0

# Get per-SNP scores
python get_snp_annots.py --fweight /storage/pandaman/project/AgentBind-GM12878-DanQ-unfixed-rnn-trans/storage/AgentBind-GM12878-DanQ/tmp/${factor}+GM12878/seqs_one_hot_c/vis-weights-total/weight.txt --resultdir /storage/mgymrek/agent-bind/singletons --TFname ${factor}
# Get regions
cat /storage/mgymrek/agent-bind/singletons/${factor}/scores.tab | \
    grep -v chrom | sed 's/chr//' | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $0}' | sort -k1,1 -k2,2n -T /storage/mgymrek/agent-bind/singletons/${factor}/ > \
    /storage/mgymrek/agent-bind/singletons/${factor}/regions.bed
# Combine with allele scores. All these commands needed to deal with off by one errors...
echo "chrom,pos,ref,alt,freq,raw.score,snr.score,core" | sed 's/,/\t/g' > /storage/mgymrek/agent-bind/singletons/${factor}/factor_singletons.tab
tabix -R /storage/mgymrek/agent-bind/singletons/${factor}/regions.bed \
    /storage/mgymrek/agent-bind/singletons/1000Genomes/1KG_afreqs.vcf.gz | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $0}' | \
    intersectBed -a stdin -b /storage/mgymrek/agent-bind/singletons/${factor}/regions.bed -wa -wb | cut -f 1-3 --complement | \
    cut -f 1,2,4-6,12,13,14,15 | \
    datamash -g 1,2,3,4,5 max 6 max 7 max 8 max 9 >> /storage/mgymrek/agent-bind/singletons/${factor}/factor_singletons.tab
