#!/bin/bash
for i in $(seq 22)
do
	chr="chrm"$i
	chrms+=($chr)
done
chrms+=("chrmX")

resolution=10kb

### Designate home directory
# Supply one mapped RNA-seq file from the previous step
##################################################
# CHANGE THIS!
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
rnaseq=$homedir/rnaseq/rnaseq_ENCFF781YWT.txt
##################################################

echo "Home directory is $homedir"
echo "Current directory is `pwd`"

## Step 1 Normalize HiC contacts and map to genes
# Normalization method is "KR"
# use python3 normalize_and_map_inter.py --help to check usage and change default arguments such as resolution
#python3 normalize_and_map_inter.py $homedir/inter $rnaseq
## example for individual chromosome pair:
echo "Step 1 Normalizing HiC and mapping to genes..."
python3 normalize_and_map_inter.py $homedir/inter $rnaseq 


## Step 2 Collapse interaction frequencies by gene pairs using four summary statistics (mean, median, max and min)
quantile=0.9
# default quantile to cutoff interaction frequency is 0.9. 
# i.e. if mean (or median, max, min) interaction frequency of a gene pair is greater than the 90% quantile of all gene pairs in that chromosome pair
# then an edge is inferred between the two genes
echo "Step 2 Collapsing interactions by gene pairs..."
Rscript genepairs_collapse_inter.R $homedir/inter "$quantile" "$resolution"
# example for individual chromosome pair: Supply the chromosome pair as the third argument, with an underscore connecting them
#Rscript genepairs_collapse_inter.R $homedir/inter "$quantile" "$resolution" "1_2"


## Step 3 (optional) Network info
#Read HiC gene pairs data into edgelists and output graph related information such as edge count and degree distribution
#This is NOT required for the PhiMRF analysis
method="all"
# method is a variable connected to the summary statistics used. In this step, all are considered
del="FALSE"
# del is a variable connected to whether to delete the linear neighboring gene pairs in the HiC network
echo "Step 3 Generating Network info..."
for chrm1 in "${chrms[@]}"
do
	for chrm2 in "${chrms[@]}"
	do
		key="$chrm1"_"$chrm2"_"$quantile"_all_TRUE_"$del"
		Rscript get_mat_info_inter.R "$key" $homedir/inter  $homedir/rnaseq/gene_names/
	done
done


# Step 4 Output adjacency matrices from the graphs
# Read HiC gene pairs data into edgelists and output adjacency matrix for each gene network
method="mean"
# method="median"
# method="min"
# method="max"
del="FALSE" 
#del=FALSE: Not deleting all linear neighbors from HiC neighbors
echo "Step 4 Generating adjacency matrices..."
for chrm1 in "${chrms[@]}"
do
	for chrm2 in "${chrms[@]}"
	do
		
		key="$chrm1"_"$chrm2"_"$quantile"_"$method"_TRUE_"$del"
		Rscript prep_mat_for_MRF_inter.R $key  $homedir/inter $homedir/rnaseq/gene_names/
	done
done

