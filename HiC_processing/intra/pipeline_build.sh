#!/bin/bash
for i in $(seq 22)
do
	chr="chrm"$i
	chrms+=($chr)
done
chrms+=("chrmX")

cellline=GM12878
resolution=10kb

## Step 0 Setup and download raw data

### Designate home directory
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/intra/
# below uses the current directory
#homedir=`dirname "$0"`
echo "Home directory is $homedir"
echo "Current directory is `pwd`"
#mkdir $homedir/raw
#mkdir $homedir/norm
#mkdir $homedir/data
#mkdir $homedir/data/y
#mkdir $homedir/info
#mkdir $homedir/genepairs
#mkdir $homedir/linear

### Download
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -P $homedir/raw
# You can also download manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
# Save it in the raw folder

### Unzip 
# only unzip the desired resolution and filter
filter=MAPQGE30
# Default read filter is MAPQ>=30, can be changed to 
# filter=MAPQG0
# See Supplemental (Extended Experimental Procedure) of Rao et al 2014 (PMID: 25497547) for details
tar -xvzf $homedir/raw/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -C $homedir/raw/ --strip-components 4 --wildcards GM12878_combined/"$resolution"_resolution_intrachromosomal/chr*/"$filter"/chr*_"$resolution".* 


## Step 1 Normalize HiC contacts and map to genes
# use python3 normalize_and_map_intra.py --help to check usage and change default arguments such as resolution
python3 normalize_and_map_intra.py $homedir ../../rnaseq_processing/ensembl_map_coding.txt


## Step 2 Collapse interaction frequencies by gene pairs using four summary statistics (mean, median, max and min)
quantile=0.9
# default quantile to cutoff interaction frequency is 0.9. 
# i.e. if mean (or median, max, min) interaction frequency of a gene pair is greater than the 90% quantile of all gene pairs in that chromosome
# then an edge is inferred between the two genes
Rscript genepairs_collapse_intra.R $homedir "$quantile" "$resolution"
# This command takes an optional fourth argument for a single chromosome, default is all chromosomes, e.g.
#Rscript genepairs_collapse_intra.R $homedir "$quantile" "$resolution" '19'




## Step 3 (optional) Network info
#Read HiC gene pairs data into edgelists and output graph related information such as edge count and degree distribution
#This is not required for the PhiMRF analysis
method="all"
# method is a variable connected to the summary statistics used. In this step, all are considered
del="FALSE"
# del is a variable connected to whether to delete the linear neighboring gene pairs in the HiC network
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quant"_all_TRUE_"$del"
	Rscript get_mat_info_intra.R "$key" $homedir  ../../rnaseq_processing/gene_names/
done


# Step 4 Output adjacency matrices from the graphs
# Read HiC gene pairs data into edgelists and output adjacency matrix for each gene network
method="mean"
# method="median"
# method="min"
# method="max"
del="FALSE" 
#del=FALSE: Not deleting all linear neighbors from HiC neighbors
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quant"_"$method"_TRUE_"$del"
	Rscript prep_mat_for_MRF_intra.R $key  $homedir ../../rnaseq_processing/gene_names 
done

