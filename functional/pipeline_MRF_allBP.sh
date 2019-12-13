#!/bin/bash
MRFfolder=/home/nzhou/hic/IMR90/work/MRF_HIC_GE/real_functional/
scriptfolder=/home/nzhou/hic/rao2014/script/goa/
datafolder=/home/nzhou/hic/rao2014/IMR90_10kb/allBP/
prop="False"
N=20
arr=()
while IFS= read -r line || [[ "$line" ]]; 
do 
	arr+=("$line")
done < $datafolder/Top"$N"_GO_terms_"$prop".txt

for GO in "${arr[@]}"
do
	key="${GO}_${prop}"
	echo $key
	Rscript $scriptfolder/save_generic_data.R  $key >$MRFfolder/out/out_"$key".info
done


