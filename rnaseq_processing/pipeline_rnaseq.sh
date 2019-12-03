#!/bin/bash
for i in $(seq 22)
do
	chr="chrm"$i
	chrms+=($chr)
done
chrms+=("chrmX")
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
## Step 1  Ensembl gene ID and coordinates 
# (Optional, only carry out this step if you are using a non-default Ensembl release)
# Default release: 90
### Download Ensembl gene ID and coordinates mapping
# To download another release, go to https://www.ensembl.org/biomart/martview/
# Select "Attributes" on the lefthand column, in the expanded table of "GENE", select "Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)" and "Gene end (bp)"; in the expanded table of "External References", select "NCBI gene ID".
# Save the downloaded file as  ../../rnaseq_processing/ensembl_map_coding.txt 


## Step 2  Generate new linear neighborhoods according to Ensembl gene coordinates
# Default provided. 
# This has nothing to do with cellline, only carry out this step if you are using a non-default Ensembl release
cutoff=10000 
# cutoff is the cutoff for linear distance
# i.e. if two genes are more than 10000 bps apart then they are not linear neighbors
python3 genepairs_collapse_linear.py  "$homedir/intra/linear" ./ensembl_map_coding.txt  "$cutoff"
#@TODO output_adjacency_linear.R

## Step 3 Write to the genename folder
# Default folder is located at ../../rnaseq_processing/gene_names/
#@TODO


## Step 4 Process RNA-seq quantification file
# The RNA-seq quantification files are downloaded from ENCODE
# Make sure the cell line of your RNAseq file agrees with the cellline for HiC file
python3 prep_rnaseq.py  ./ensembl_map_coding.txt  "$homedir/rnaseq/ENCFF*.tsv"  $homedir/rnaseq/


# Step 5 Output y 
# y needs to be in an order that correpsonds to the adjacency matrix from Step 5
# More than one processed RNA-seq file can be supplied at the same time using wildcards
# Each RNA-seq file will be read in as a column in a matirx of rnaseq data
# NO quotes around the last argument if you are using wildcard
for chrm in "${chrms[@]}"
do
	Rscript prep_y_for_MRF_intra.R "$chrm" $homedir/intra/  ./gene_names/  $homedir/rnaseq/rnaseq_*.txt
done





