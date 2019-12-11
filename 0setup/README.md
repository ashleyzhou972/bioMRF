#!/bin/bash

cellline=GM12878
resolution=10kb

## Step 0 Setup and download raw data

### Designate home directory

###################################################
#CHANGE THIS!
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
####################################################

echo "Home directory is $homedir"
echo "Current directory is `pwd`"

# Create folders for RNA-seq data
mkdir $homedir/rnaseq
mkdir $homedir/rnaseq/gene_names

# Create folders for intra-chromosomal and inter-chromosomal
for int in 'intra', 'inter'
do
	mkdir $homedir/$int
	mkdir $homedir/$int/raw
	mkdir $homedir/$int/norm
	mkdir $homedir/$int/data
	mkdir $homedir/$int/data/y
	mkdir $homedir/$int/info
	mkdir $homedir/$int/genepairs
	mkdir $homedir/$int/linear
done


### Download intra
echo "Step 0 Downloading intra-chromosomal HiC data..."
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -P $homedir/intra/raw
# You can also download manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
# Save it in the raw folder

### Unzip intra
# only unzip the desired resolution and filter
filter=MAPQGE30
norm=KR
# Default read filter is MAPQ>=30, can be changed to 
# filter=MAPQG0
# See Supplemental (Extended Experimental Procedure) of Rao et al 2014 (PMID: 25497547) for details
echo "Step 0 Extracting intra-chromosomal HiC data"
tar -xvzf $homedir/intra/raw/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -C $homedir/intra/raw/ --strip-components 4 --wildcards GM12878_combined/"$resolution"_resolution_intrachromosomal/chr*/"$filter"/*.{RAWobserved,"$norm"norm}

### Download inter
echo "Step 0 Downloading inter-chromosomal HiC data..."
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_combined_interchromosomal_contact_matrices.tar.gz -P $homedir/inter/raw
# You can also download manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
# Save it in the raw folder

### Unzip inter
# only unzip the desired resolution and filter
filter=MAPQGE30
norm=KR
# Default read filter is MAPQ>=30, can be changed to 
# filter=MAPQG0
# See Supplemental (Extended Experimental Procedure) of Rao et al 2014 (PMID: 25497547) for details
echo "Step 0 Extracting inter-chromosomal HiC data..."
tar -xvzf $homedir/inter/raw/GSE63525_"$cellline"_combined_interchromosomal_contact_matrices.tar.gz -C $homedir/inter/raw/ --strip-components 4 --wildcards GM12878_combined_interchromosomal/"$resolution"_resolution_interchromosomal/chr*_chr*/"$filter"/*.{RAWobserved,INTER"$norm"norm}


## Step 1  Download Ensembl gene ID and coordinates 
# (Optional, only carry out this step if you are using a non-default Ensembl release)
# Default release: 90
### Download Ensembl gene ID and coordinates mapping
# To download another release, go to https://www.ensembl.org/biomart/martview/
# Select "Attributes" on the lefthand column, in the expanded table of "GENE", select "Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)" and "Gene end (bp)"; in the expanded table of "External References", select "NCBI gene ID".
# Save the downloaded file as  ../../1rnaseq_processing/ensembl_map_coding.txt 


## Step 2 Download RNAseq data and map to Ensembl gene ID
# The RNA-seq quantification files are downloaded from ENCODE.
# Make sure the cell line of your RNAseq file agrees with the cellline for HiC file
# If you wish to use non-ENCODE RNAseq quantification file, make sure the file follows RSEM's output format:
# https://www.encodeproject.org/documents/0c78ea4b-9392-421b-a6f3-6c858b6002aa/@@download/attachment/RSEM_Documentation.pdf
# NOTE that the first column should be gene_id and second column is transcript_id, contrary to the above document.
# The quantification metric used in this project is the 8th column: posterior_mean_count
# If you wish to use another COUNT metric, please specify the column index (1-based) in the next pipeline for rnaseq.

# Example for GM12878 (two isogenic replicates)
wget https://www.encodeproject.org/files/ENCFF781YWT/@@download/ENCFF781YWT.tsv -P $homedir/rnaseq/
wget https://www.encodeproject.org/files/ENCFF680ZFZ/@@download/ENCFF680ZFZ.tsv -P $homedir/rnaseq/


