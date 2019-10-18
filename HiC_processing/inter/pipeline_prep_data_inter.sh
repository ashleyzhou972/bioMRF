#!/bin/bash
#Step 1 Collapse normalized HiC counts to gene level by the four methods (mean, median, max and min)
#Rscript get_p_values_inter.R

#Step 2 Read these HiC data into edgelists and output graph related information such as edge count and degree
method="all"
del='FALSE'
for chrm1 in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
#for chrm1 in 'chrm1' 'chrm7' 'chrm14' 'chrm17' 'chrmX'
do 
	for chrm2 in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
	do
		for quant in '0.9'
		do		
			key="$chrm1"_"$chrm2"_"$quant"_"$method"_TRUE_"$del"
			#echo $key
			Rscript get_mat_info_inter.R $key 
		done
	done
done
#Step 3 Output adjacency matrices from the graphs
del='FALSE'
for chrm1 in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
#for chrm in 'chrm1'
#for chrm1 in 'chrm1' 'chrm7' 'chrm14' 'chrm17' 'chrmX'
do 
	for chrm2 in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
	#for chrm2 in 'chrm1' 'chrm7' 'chrm14' 'chrm17' 'chrmX'
	do
		for quant in '0.9'
		do		
			for method in 'mean' 'median' 'max' 'min'
			do
				key="$chrm1"_"$chrm2"_"$quant"_"$method"_TRUE_"$del"
				#echo $key
				Rscript prep_mat_for_MRF_inter.R $key 

			done
		done
	done
done

#Step 4 Output y (RNAseq) data for each chromosome according to the order of genes in the adjacency matrices (One y per chromosome regardless of HiC network parameters)
for chrm1 in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
#for chrm in 'chrm1'
#for chrm1 in 'chrm1' 'chrm7' 'chrm14' 'chrm17' 'chrmX'
do 
	for chrm2 in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
	#for chrm2 in 'chrm1' 'chrm7' 'chrm14' 'chrm17' 'chrmX'
	do
		Rscript prep_y_for_MRF_inter.R "$chrm1" "$chrm2"
	done
done



