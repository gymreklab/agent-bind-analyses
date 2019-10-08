#!/bin/bash

TRAIT=$1

SUMM=/storage/mgymrek/agent-bind/ldsc/summstats/${TRAIT}.sumstats.gz
#BASELINE=/storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_baseline/baseline.
BASELINE=/storage/mgymrek/agent-bind/ldsc/resources/1000G_EUR_Phase3_baseline_v2.0.2/baselineLD.
WEIGHTS=/storage/mgymrek/agent-bind/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=/storage/mgymrek/agent-bind/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC.

# LD scores for each tool
ABLOCUS=/storage/mgymrek/agent-bind/ldsc/annotations/agentbind-locus-quant/AgentBindLocusQuant.
#/storage/mgymrek/agent-bind/ldsc/annotations/agentbind-locus/AgentBindLocus.
IMPACT_FOXP3=/storage/mgymrek/workspace/IMPACT/sLDSC/IMPACTannotations/Foxp3_GenomewideTrack_IMPACT.
IMPACT_GATA3=/storage/mgymrek/workspace/IMPACT/sLDSC/IMPACTannotations/Gata3_GenomewideTrack_IMPACT.
IMPACT_STAT3=/storage/mgymrek/workspace/IMPACT/sLDSC/IMPACTannotations/Stat3_GenomewideTrack_IMPACT.

python /storage/resources/source/ldsc/ldsc.py \
    --h2 $SUMM \
    --ref-ld-chr $BASELINE,${IMPACT_FOXP3},${IMPACT_GATA3},${IMPACT_STAT3},$ABLOCUS \
    --w-ld-chr $WEIGHTS \
    --overlap-annot \
    --frqfile-chr $FRQ \
    --print-coefficients --print-delete-vals \
    --out /storage/mgymrek/agent-bind/ldsc/results/${TRAIT}
