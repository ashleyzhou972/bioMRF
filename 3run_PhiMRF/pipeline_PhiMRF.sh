#!/bin/bash

for i in $(seq 22)
do
	chr="chrm"$i
	chrms+=($chr)
done
chrms+=("chrmX")

#############################################
# CHANGE THIS!
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
quantile=0.9
method="mean"
############################################


mkdir $homedir/intra/results/
mkdir $homedir/inter/results


for chrmA in "${chrms[@]}"
#for chrmA in "chrm2"
do
	echo "Beginning intra" "$chrmA"
	Rscript run_PhiMRF_intra.R $homedir/intra "$chrmA" "$quantile" "$method"
	echo "Beginning linear" "$chrmA"
	Rscript run_PhiMRF_linear.R $homedir/intra "$chrmA" 
	for chrmB in "${chrms[@]}"
	#for chrmB in "chrm3"
	do
		echo "Beginning inter" "$chrmA" "$chrmB"
		Rscript run_PhiMRF_inter.R $homedir/inter "$chrmA" "$chrmB" "$quantile" "$method"
	done
done
