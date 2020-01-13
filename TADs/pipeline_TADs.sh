#!/bin/bash

#######################################
#CHANGE THIS!
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
rnaseq=$homedir/rnaseq/rnaseq_ENCFF781YWT.txt
cellline=GM12878_primary+replicate
#######################################

# Check the download file name for other cell lines
# For example
#cellline=IMR90

### Step 0 set up
mkdir $homedir/TADs
mkdir $homedir/TADs/TADgenes
mkdir $homedir/TADs/TADs_all_data
mkdir $homedir/TADs/TADs_all_data/y
mkdir $homedir/TADs/TADs_non_data
mkdir $homedir/TADs/TADs_non_data/y
mkdir $homedir/TADs/TADs_intra_data
mkdir $homedir/TADs/TADs_intra_data/y
mkdir $homedir/TADs/TADs_inter_data
mkdir $homedir/TADs/TADs_inter_data/y

### Step 1 download domain list
echo "Step 1 Download TAD list..."
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_Arrowhead_domainlist.txt.gz -P $homedir/TADs
gunzip $homedir/TADs/GSE63525_"$cellline"_Arrowhead_domainlist.txt.gz


### Step 2 Map genes to TADs
echo "Step 2 Map genes to TADs..."
python3 map_ensembl_genes_tad.py $homedir/TADs/GSE63525_"$cellline"_Arrowhead_domainlist.txt $rnaseq $homedir/TADs/TADgenes
# optional: specify one chromosome
#python3 map_ensembl_genes_tad.py $homedir/TADs/GSE63525_"$cellline"_Arrowhead_domainlist.txt $rnaseq $homedir/TADs/TADgenes --chr 13


### Step 3 Save gene networks for TAD genes
# Four categories:
# - TADs_all_data: all TAD genes, all edges
# - TADs_non_data: all nonTAD genes, all edges
# - TADs_intra_data: all TAD genes, only intra-TAD edges (edges that connect two genes in the same TAD)
# - TADs_inter_data: all TAD genes, only inter-TAD edges (edges that connect two genes not in the same TAD, still in the same chromosome)

# Note that the following three step requires the processed adjacency matrices from the full intra-chromosomal networks.
# Please check the $homedir/intra/data folder to make sure it is non-empty and has by-chromosome adjacency matrices and y's.

#### Step 3.1
echo "Step 3.1 save TAD gene network..."
Rscript save_TAD_data_allTAD.R $homedir/TADs/TADgenes $homedir/intra/data/ $homedir/TADs/TADs_all_data/
#### Step 3.2
echo "Step 3.2 save non-TAD gene network..."
Rscript save_TAD_data_nonTAD.R $homedir/TADs/TADgenes $homedir/intra/data/ $homedir/TADs/TADs_non_data/
#### Step 3.3
echo "Step 3.3 save TAD gene network with intraTAD edges..."
Rscript save_TAD_data_intraTAD.R $homedir/TADs/TADgenes $homedir/intra/data/ $homedir/TADs/TADs_intra_data/

# The following step requires the TADs_all_data and TADs_intra_data folders to be filled, as described above (Step 3.1 and 3.3)
# Note the change in the second argument 
#### Step 3.4
echo "Step 3.3 save TAD gene network with interTAD edges..."
Rscript save_TAD_data_interTAD.R $homedir/TADs/TADgenes $homedir/TADs/TADs_all_data/ $homedir/TADs/TADs_inter_data/



