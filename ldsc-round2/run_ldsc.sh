#!/bin/bash

SUMM=/storage/mgymrek/agent-bind/ldsc/summstats/RA_GWASmeta_European_v2.sumstats.gz
BASELINE=/storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_baseline/baseline.
WEIGHTS=/storage/mgymrek/agent-bind/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=/storage/mgymrek/agent-bind/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC.

python /storage/resources/source/ldsc/ldsc.py \
    --h2 $SUMM \
    --ref-ld-chr $BASELINE \
    --w-ld-chr $WEIGHTS \
    --overlap-annot \
    --frqfile-chr $FRQ \
    --out /storage/mgymrek/agent-bind/ldsc/results/RA_GWASmeta_European_v2
