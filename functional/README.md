# Process functional data

## **Note!** The [biopython](https://biopython.org) package is required!

All of the following steps are in [`pipeline_functional.sh`](pipeline_functional.sh) 

Please read through each step and customize, before running

```
bash pipeline_functional.sh
```

### Designate home dirctory

The same directory as the [set up](../0setup/) step.

Example:
```
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
```

### Supply one processed RNA-seq file from the previous step
Example:
```
rnaseq=$homedir/rnaseq/rnaseq_ENCFF781YWT.txt
```

### Supply the obo file name
Example:
```
oboname=go.obo.gz
```

### Step 0 set up

```
mkdir $homedir/functional
mkdir $homedir/functional/data
mkdir $homedir/functional/data/y
mkdir $homedir/functional/annotations
mkdir $homedir/functional/info
```

### Step 1 download Gene Ontology Annotations
 1. Download obo file
See [Gene Ontology Consortium](ftp://ftp.geneontology.org/go/www/GO.format.obo-1_2.shtml) for details of the obo format

To download the release used in the paper (20190701):
```
wget ftp://ftp.geneontology.org/go/ontology-archive//gene_ontology_edit.obo.2019-07-01.gz -O $homedir/functional/annotations/$oboname
```
To download the most recent release:
```
wget http://purl.obolibrary.org/obo/go.obo -O $homedir/functional/annotations/$oboname
```

 2. Unzip go.obo
```
gunzip $homedir/functional/annotations/$oboname
```

 3. Download GOA

To download the release used in the paper (20190916):
```
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old//HUMAN/goa_human.gaf.194.gz -O $homedir/functional/annotations/goa_human.gaf.gz
```

To download the most recent release:
```
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa//HUMAN/goa_human.gaf.gz -P $homedir/functional/annotations/goa_human.gaf.gz
```

 4. Unzip gaf file
```
gunzip $homedir/functional/annotations/goa_human.gaf.gz
```

 5. Download mapping between uniprot and Ensembl

Note you can change which Ensembl release to use
```
releasenum=90
wget ftp://ftp.ensembl.org/pub/release-"$releasenum"/tsv/homo_sapiens//Homo_sapiens.GRCh38."$releasenum".uniprot.tsv.gz -O $homedir/functional/annotations/Ensembl_uniprot.tsv.gz
```

### Step 2 Get the most annotated proteins (genes) in Biological Process Ontology (Optional!) 
This step is optional, the list of top 20 GO terms is already provided as `Top20_GO_terms_counts_False` (unpropagated) and `Top20_GO_terms_counts_True` (propagated) in this directory
List of Ensembl gene names annotated with each GO term is also already provided as the `ensmebl_list` folder in this directory
Below code outputs the same files as above, but in the $homedir directory

Check usage (`python3 get_bpo_term_with_most_annotations.py --help`) to see options to 

 - Output the top N GO terms
 - propagate the annotations through the GO hierarchy (including MFO terms)
 - exclude direct children of root terms

```
prop="False"
N=20
python3 get_bpo_term_with_most_annotations.py --no-propagate $N $homedir/functional/annotations/  $homedir/rnaseq/gene_names/
arr=()

while IFS= read -r line || [[ "$line" ]]; 
do 
	term="$(cut -f1 <<<"$line")"
	echo $term
	arr+=("$term")
done < ./Top"$N"_GO_terms_counts_"$prop"

```

### Step 3 Aggregate all intra- and inter-chromosomal adjacency matrices

This could take a while...
```
Rscript combine_nets.R $homedir
```

### Step 4 Output data for subsets of genes for each GO term
```
for GO in "${arr[@]}"
do
	key="${GO}_${prop}"
	echo $key
	Rscript save_generic_data.R  $key >$homedir/functional/info/out_"$key".info
done
```
