#!/bin/bash

EXPERIMENT=ICL103
OMICS=Proteins
VARIABLE=VN1203
STRAIN="VN1203 NS1 627E Mock"
BASELINE=NS1
TIME="0h 3h 7h 12h 18h 24h"


for CONSTANT in $TIME
do
    Rscript ProgenyDorotheaInfluenza.R  --variable $VARIABLE \
	    --baseline $BASELINE \
	    --constant $CONSTANT \
	    --organism Human \
	    --outdir  IntermediateFiles \
	    Dorothea_${EXPERIMENT}_${OMICS}_${VARIABLE}_vs_${BASELINE}_${CONSTANT}r.csv \
	    Progeny_${EXPERIMENT}_${OMICS}_${VARIABLE}_vs_${BASELINE}_${CONSTANT}r.csv 
done
