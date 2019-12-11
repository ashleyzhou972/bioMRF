#!/bin/bash
for i in $(seq 22)
do
	chr="chrm"$i
	chrms+=($chr)
done
chrms+=("chrmX")
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/

## Step 0 
echo "Please make sure you have downloaded/prepared your RNA-seq quantification file."
# Please refer to the instructions in ../0setup/pipeline_setup, step 2.

## Step 1 Process RNA-seq quantification file
# The RNA-seq quantification files are downloaded from ENCODE
echo "Step 1 Process RNA-seq quantification file..."
############################################################
# CHANGE BELOW!!
# replace the second argument with your RAW RNA-seq filename
python3 prep_rnaseq.py  ./ensembl_map_coding.txt  "$homedir/rnaseq/ENCFF*.tsv"  $homedir/rnaseq/
# designate the path for the PROCESSED RNA-seq file
# only one file is enough
# naming convention is "rnaseq_[Your raw RNA-seq BASE filename].txt"
# CHANGE BELOW!!
rnaseq=$homedir/rnaseq/rnaseq_ENCFF781YWT.txt
############################################################


## Step 2 Write to the genename folder
# Our model should run on all genes for which RNA-seq data is available.
# Therefore, once we've mapped all the RNA-seq data we write all the names of the genes by chromosome to a separate folder
# To inform subsequent mapping of HiC data
echo "Step 2 Write gene names..."
Rscript save_gene_names.R $rnaseq  $homedir/rnaseq/gene_names/


# Step 3 Output y 
# y needs to be in an order that correpsonds to the adjacency matrix from Step 5
# More than one processed RNA-seq file can be supplied at the same time using wildcards
# Each RNA-seq file will be read in as a column in a matirx of rnaseq data
# NO quotes around the last argument if you are using wildcard
echo "Step 3 Ouput y for PhiMRF (both intra and inter)..."
for chrmA in "${chrms[@]}"
do
	Rscript prep_y_for_MRF_intra.R "$chrmA" $homedir/intra/  $homedir/rnaseq/gene_names/  $homedir/rnaseq/rnaseq_*.txt
	for chrmB in "${chrms[@]}"
	do
		Rscript prep_y_for_MRF_inter.R "$chrmA" "$chrmB" $homedir/inter/  $homedir/rnaseq/gene_names/  $homedir/rnaseq/rnaseq_*.txt
	done
done





