#Step 1 Collapse normalized HiC counts to gene level by the four methods (mean, median, max and min)
#Rscript get_p_values_intra.R

#Step 2 Read these HiC data into edgelists and output graph related informationsuch as edge count and degree
for chrm in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
#for chrm in 'chrm9'
do
	for quant in '0.9'
	do		
		method="all"
		for del in 'TRUE' 'FALSE'
		do
			key="$chrm"_"$quant"_"$method"_TRUE_"$del"
			#echo $key
			Rscript get_mat_info_intra.R $key 
		done
	done
done
#Step 3 Output adjacency matrices from the graphs
for chrm in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
#for chrm in 'chrm9'
do
	for quant in '0.9'
	do		
		for method in 'mean' 'median' 'max' 'min'
		do
			for del in 'TRUE' 'FALSE'
			do
				key="$chrm"_"$quant"_"$method"_TRUE_"$del"
				#echo $key
				Rscript prep_mat_for_MRF_intra.R $key 
			done

		done
	done
done

#Step 4 Output y (RNAseq) data for each chromosome according to the order of genes in the adjacency matrices (One y per chromosome regardless of HiC network parameters)
for chrm in 'chrm1' 'chrm2' 'chrm3' 'chrm4' 'chrm5' 'chrm6' 'chrm7' 'chrm8' 'chrm9' 'chrm10' 'chrm11' 'chrm12' 'chrm13' 'chrm14' 'chrm15' 'chrm17' 'chrm16' 'chrm18' 'chrm19' 'chrm20' 'chrm21' 'chrm22' 'chrmX'
#for chrm in 'chrm9'
do
	Rscript prep_y_for_MRF_intra.R "$chrm"
done
