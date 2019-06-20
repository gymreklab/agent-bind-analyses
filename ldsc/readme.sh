#!/bin/bash

source activate ldsc
source params.sh
CHROMS=$(seq 1 22)

#./get_abscores_all.sh

./get_impact.sh 

# Run LDSC using annotations above

python /storage/cynthiawu/soft/ldsc/ldsc.py \
    --h2 ${SUMSTATS} \
    --ref-ld-chr ${OUTDIR}/annot/IMPACT.,${BASELINE} \
    --w-ld-chr ${WEIGHTS} \
    --overlap-annot \
    --frqfile-chr ${FREQ} \
    --out test.IMPACT
    
for TF in Gata3 Foxp3 Tbet Stat3
do
    python /storage/cynthiawu/soft/ldsc/ldsc.py \
	--h2 ${SUMSTATS} \
	--ref-ld-chr ${OUTDIR}/annot/${TF}.AB.,${BASELINE} \
	--w-ld-chr ${WEIGHTS} \
	--overlap-annot \
	--frqfile-chr ${FREQ} \
	--out test.${TF}
done
