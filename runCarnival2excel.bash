#!/bin/bash

EXPERIMENT=ICL103
OMICS=Proteins
VARIABLE=VN1203
STRAIN="VN1203 NS1 627E Mock"
BASELINE=NS1
TIME="0h 3h 7h 12h 18h 24h"
OMNIPATHDB=OmnipathSignedDirectedInteractions.csv
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


for BASELINE in Mock NS1 
do
    if [ $BASELINE = "NS1" ]
    then
	PERTURBATION=HostVirus
    else
	PERTURBATION=NoInput
    fi
    
    for CONSTANT in $TIME
    do
	INPUTDIR=ResultsCARNIVAL_${VARIABLE}vs${BASELINE}_for_${CONSTANT}r
	OUTPUTDIR=ResultsCARNIVAL_${VARIABLE}vs${BASELINE}_for_${CONSTANT}r
	echo $INPUTDIR
	if [ -r ${INPUTDIR}/network_solution.dot ]
	   then
	       python dot2excel.py \
		      --perturbation_file $PERTURBATION \
		      --variable $VARIABLE \
		      --baseline $BASELINE \
		      --constant $CONSTANT \
		      --inputdir $INPUTDIR \
		      --outputdir $OUTPUTDIR \
		      --omnipathdb $OMNIPATHDB
	       
	       dot -Tpdf \
		   ${INPUTDIR}/network_solution.dot \
		   -o ${OUTPUTDIR}/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}.pdf
	fi
	
    done
done
