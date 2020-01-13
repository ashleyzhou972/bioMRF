# Process RNA-seq data

All of the following steps are in [`pipeline_rnaseq.sh`](pipeline_rnaseq.sh). 

Please read through each step and customize, before running

`bash pipeline_rnaseq.sh`

### Designate home dirctory

The same directory as the [set up](../0setup/) step.

Example:
```
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
```

## Process RNA-seq quantification file

Make sure you have downloaded/prepared your RNA-seq files, according to the [set up](../0setup/) step.

1. Replace the second argument below with your ***raw*** RNA-seq file path, wildcards accepted.
```
python3 prep_rnaseq.py  ./ensembl_map_coding.txt  "$homedir/rnaseq/ENCFF*.tsv"  $homedir/rnaseq/
```
2. Supply the path for the ***processed*** RNA-seq file below. Only one of the files is enough.

Naming convention for processed files is "rnaseq_[Your raw RNA-seq BASE filename].txt".

Example:
```
rnaseq=$homedir/rnaseq/rnaseq_ENCFF781YWT.txt
```


### Write to the genename folder

PhiMRF should use all genes for which RNA-seq data is available. Therefore, once we've mapped all the RNA-seq data, we write all the names of the genes by chromosome to a separate folder, to inform subsequent mapping of HiC data. 
```
Rscript save_gene_names.R $rnaseq  $homedir/rnaseq/gene_names/
```


### Output RNA-seq count 

This step outputs the RNA-seq count in an R matrix format, each RNA-seq file (replicate) will be one column in the matrix.

More than one processed RNA-seq file can be supplied at the same time using wildcards. (NO quotes around the last argument if you are using wildcard)

Although the gene expression data are the same, we output two set of files, for intra- and inter-chromosomal datasets.

```
for chrmA in "${chrms[@]}"
do
	Rscript prep_y_for_MRF_intra.R "$chrmA" $homedir/intra/  $homedir/rnaseq/gene_names/  $homedir/rnaseq/rnaseq_*.txt
	for chrmB in "${chrms[@]}"
	do
		Rscript prep_y_for_MRF_inter.R "$chrmA" "$chrmB" $homedir/inter/  $homedir/rnaseq/gene_names/  $homedir/rnaseq/rnaseq_*.txt
	done
done
```

### Go to the [next step](../2hic_processing)
