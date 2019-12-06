#!/bin/bash

for i in $(seq 22)
do
	chr="chrm"$i
	chrms+=($chr)
done
chrms+=("chrmX")

homedir=/home/nzhou/hic/rao2014/GM12878_10kb/intra/
mkdir $homedir/results/

#for chrm in "${chrms[@]}"
for chrm in "chrm1"
do
	Rscript example_PhiMRF.R $homedir "$chrm" 
done
