#!/bin/bash

#######################################
#CHANGE THIS!
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
rnaseq=$homedir/rnaseq/rnaseq_ENCFF781YWT.txt
oboname=go.obo.gz
#######################################

### Step 0 set up
mkdir $homedir/functional
mkdir $homedir/functional/data
mkdir $homedir/functional/data/y
mkdir $homedir/functional/annotations
mkdir $homedir/functional/info

### Step 1 download Gene Ontology Annotations
echo "Step 1.1 Download go.obo..."
## See ftp://ftp.geneontology.org/go/www/GO.format.obo-1_2.shtml for details of the obo format
## Release from 20190701:
#wget ftp://ftp.geneontology.org/go/ontology-archive//gene_ontology_edit.obo.2019-07-01.gz -O $homedir/functional/annotations/$oboname
## Below downloads the most release
#wget http://purl.obolibrary.org/obo/go.obo -O $homedir/functional/annotations/$oboname

echo "Unzip go.obo..."
#gunzip $homedir/functional/annotations/$oboname

echo "Step 1.2 Download GOA..."
# Release from 20190916
#wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old//HUMAN/goa_human.gaf.194.gz -O $homedir/functional/annotations/goa_human.gaf.gz

# Below downloads the most recent release
#wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa//HUMAN/goa_human.gaf.gz -P $homedir/functional/annotations/goa_human.gaf.gz

echo "Unzip gaf file..."
#gunzip $homedir/functional/annotations/goa_human.gaf.gz

echo "Step 1.3 Download mapping between uniprot and Ensembl..."
# Note you can change which Ensembl release to use
releasenum=90
#wget ftp://ftp.ensembl.org/pub/release-"$releasenum"/tsv/homo_sapiens//Homo_sapiens.GRCh38."$releasenum".uniprot.tsv.gz -O $homedir/functional/annotations/Ensembl_uniprot.tsv.gz


### Step 2 Get the most annotated proteins (genes) in Biological Process Ontology (Optional!) 
## This step is optional, the list of top 20 GO terms is already provided as "Top20_GO_terms_counts_False" (unpropagated) and "Top20_GO_terms_counts_True" (propagated) in this directory
## List of Ensembl gene names annotated with each GO term is also already provided as the "ensmebl_list" folder in this directory
## Below code outputs the same files as above, but in the $homedir directory

## Check usage to see options to 
# - Output the top N GO terms
# - propagate the annotations through the GO hierarchy (including MFO terms)
# - exclude direct children of root terms

echo "Step 2 Get BPO terms with the most annotations..."
prop="False"
N=20
#python3 get_bpo_term_with_most_annotations.py --no-propagate $N $homedir/functional/annotations/  $homedir/rnaseq/gene_names/
arr=()
while IFS= read -r line || [[ "$line" ]]; 
do 
	term="$(cut -f1 <<<"$line")"
	echo $term
	arr+=("$term")
done < ./Top"$N"_GO_terms_counts_"$prop"

### Step 3 Aggregate all intra- and inter-chromosomal adjacency matrices (This could take a while...)
echo "Step3 3 Aggregate all intra- and inter- adjacency matrices..."
#Rscript combine_nets.R $homedir

### Step 4 Output data for subsets of genes for each GO term
echo "Step 4 Save data for PhiMRF for each GO term..."

#use the array saved in Step 2
for GO in "${arr[@]}"
do
	key="${GO}_${prop}"
	echo $key
	Rscript save_generic_data.R  $key $homedir
done
