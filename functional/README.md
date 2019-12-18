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

### Step 0 set up

```
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
```

### Step 1 download domain list

Here we used the TADs called by the Arrowhead caller([Rao et al 2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547)), downloaded from [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) (See Pages S50-S51 in the Supplemental 1 (Extended Experimental Procedure) of [Rao et al 2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547) for how Arrowhead works).

A number of other TAD callers are available and may output different results. See [Forcatto et al 2017](https://www.ncbi.nlm.nih.gov/pubmed/28604721) for a summary of other available TAD callers.

```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_Arrowhead_domainlist.txt.gz -P $homedir/TADs
gunzip $homedir/TADs/GSE63525_"$cellline"_Arrowhead_domainlist.txt.gz
```


### Step 2 Map genes to TADs

```
python3 map_ensembl_genes_tad.py $homedir/TADs/GSE63525_"$cellline"_Arrowhead_domainlist.txt $rnaseq $homedir/TADs/TADgenes
```
optionally, you can specify only one chromosome, `python3 map_ensembl_genes_tad.py $homedir/TADs/GSE63525_"$cellline"_Arrowhead_domainlist.txt $rnaseq $homedir/TADs/TADgenes --chr 13`.


### Step 3 Save gene networks for TAD genes

Four categories:

  - TADs_all_data: all TAD genes, all edges
  - TADs_non_data: all nonTAD genes, all edges
  - TADs_intra_data: all TAD genes, only intra-TAD edges (edges that connect two genes in the same TAD)
  - TADs_inter_data: all TAD genes, only inter-TAD edges (edges that connect two genes not in the same TAD, still in the same chromosome)

Note that the following 3 step requires the processed adjacency matrices from the full intra-chromosomal networks.

Please check the `$homedir/intra/data` folder to make sure it is non-empty and has by-chromosome adjacency matrices and y's.

```
Rscript save_TAD_data_allTAD.R $homedir/TADs/TADgenes $homedir/intra/data/ $homedir/TADs/TADs_all_data/
Rscript save_TAD_data_nonTAD.R $homedir/TADs/TADgenes $homedir/intra/data/ $homedir/TADs/TADs_non_data/
Rscript save_TAD_data_intraTAD.R $homedir/TADs/TADgenes $homedir/intra/data/ $homedir/TADs/TADs_intra_data/
```

The following step requires the `$homedir/TADs/TADs_all_data` and `$homedir/TADs/TADs_intra_data` folders to be filled, as described above.

Note the change in the second argument below.
```
Rscript save_TAD_data_interTAD.R $homedir/TADs/TADgenes $homedir/TADs/TADs_all_data/ $homedir/TADs/TADs_inter_data/
```
