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
# use python3 normalize_and_map_intra.py --help to check usage and change default arguments such as resolution
echo "Step 1 Normalizing HiC and mapping to genes..."
python3 normalize_and_map_intra.py $homedir/intra $rnaseq --resolution "$resolution" --overlap 0.1 --norm 'KR


## Step 2 Collapse interaction frequencies by gene pairs using four summary statistics (mean, median, max and min)
quantile=0.9
# default quantile to cutoff interaction frequency is 0.9. 
# i.e. if mean (or median, max, min) interaction frequency of a gene pair is greater than the 90% quantile of all gene pairs in that chromosome
# then an edge is inferred between the two genes
echo "Step 2 Collapsing interactions by gene pairs..."
Rscript genepairs_collapse_intra.R $homedir/intra "$quantile" "$resolution"
# This command takes an optional fourth argument for a single chromosome, default is all chromosomes, e.g.
#Rscript genepairs_collapse_intra.R $homedir/intra "$quantile" "$resolution" '19'


## Step 3 Generate linear neighborhoods according to Ensembl gene coordinates
cutoff=10000 
# cutoff is the cutoff for linear distance
# i.e. if two genes are more than 10000 bps apart then they are not linear neighbors
echo "Step 3 Linear neighborhoods..."
python3 genepairs_collapse_linear.py  "$homedir/intra/linear" $rnaseq  "$cutoff"
Rscript prep_mat_for_MRF_linear.R $homedir/intra/linear $rnaseq $homedir/rnaseq/gene_names/ "TRUE"





## Step 4 (optional) Network info
#Read HiC gene pairs data into edgelists and output graph related information such as edge count and degree distribution
#This is not required for the PhiMRF analysis
method="all"
# method is a variable connected to the summary statistics used. In this step, all are considered
del="FALSE"
# del is a variable connected to whether to delete the linear neighboring gene pairs in the HiC network
echo "Step 4 Generating Network info..."
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quantile"_all_TRUE_"$del"
	Rscript get_mat_info_intra.R "$key" $homedir/intra  $homedir/rnaseq/gene_names/
done


# Step 5 Output adjacency matrices from the graphs
# Read HiC gene pairs data into edgelists and output adjacency matrix for each gene network
method="mean"
# method="median"
# method="min"
# method="max"
del="FALSE" 
#del=FALSE: Not deleting all linear neighbors from HiC neighbors
echo "Step 5 Generating adjacency matrices..."
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quantile"_"$method"_TRUE_"$del"
	Rscript prep_mat_for_MRF_intra.R $key  $homedir/intra $homedir/rnaseq/gene_names 
done


