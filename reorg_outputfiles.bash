VARIABLE=VN1203
STRAIN="Mock NS1"
BASELINE=Mock
TIME="0h  12h 18h 24h"
PERTURBATION=kegg_influenza_ns1
for BASELINE in $STRAIN
do

    for CONSTANT in $TIME
    do
	INPUTDIR=Results/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}
	git mv $INPUTDIR/network_solution.dot Results/Dotfiles/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}.dot
	git mv $INPUTDIR/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}.pdf Results/Networks/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}.pdf
	#git mv Results/Networks/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}.xlsx  Results/ExcelFiles/${VARIABLE}vs${BASELINE}_for_${CONSTANT}r_${PERTURBATION}.xlsx
    done
done
git status
