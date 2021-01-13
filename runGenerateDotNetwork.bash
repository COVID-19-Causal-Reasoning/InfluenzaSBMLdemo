#!/bin/bash

EXPERIMENT=ICL103
OMICS=Proteins
VARIABLE=VN1203
STRAIN="VN1203 NS1 627E Mock"
BASELINE=Mock
TIME="0h 3h 7h 12h 18h 24h"
PERTURBATION=kegg_influenza_ns1


# for CONSTANT in $TIME
# do
#     mkdir -p ResultsCARNIVAL_${VARIABLE}vs${BASELINE}_for_${CONSTANT}r 
#     Rscript runCarnival.R  --variable $VARIABLE \
# 	    --baseline $BASELINE \
# 	    --constant $CONSTANT \
# 	    --threads 8 \
# 	    --outdir  ResultsCARNIVAL_${VARIABLE}vs${BASELINE}_for_${CONSTANT}r \
# 	    --perturbation_file $PERTURBATION \
# 	    IntermediateFiles/${VARIABLE}vs${BASELINE}_for_${CONSTANT}_pathways_activity_score_inputCarnival.csv \
# 	    IntermediateFiles/${VARIABLE}vs${BASELINE}_for_${CONSTANT}_tf_activities_stat_inputCarnival.csv 
# done


for BASELINE in Mock 
do
    for CONSTANT in $TIME
    do
	OUTDIR=ResultsCARNIVAL${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}
	dot -Tpdf  ${OUTDIR}/network_solution.dot -o ${OUTDIR}/${OUTDIR}.pdf
    done
done
