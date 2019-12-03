# bioMRF: Exploring the Spatial Dependency of Gene Expression
This repository contains 
- scripts to *organize and clean* HiC and RNAseq data
- scripts to run the [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) model on RNAseq and HiC data
- The source code of the R package [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) is available as a separate repo.


The pipeline can be separated into three parts, processing HiC data, processing RNAseq data and running PhiMRF.
The three pipelines are explained below, as well as in their separate `pipeline_*.sh` files (In HiC processing, the pipeline is further separated into intra- and inter-chromosoma).

If you would like run the pipelines from scratch, simply run `bash pipeline_*.sh` in their respective folders.

However, some results in those pipelines are optional or already provided in this repo, so it is recommended that you read through the pipelines, understand each step, and comment out the steps that are unnecessary for your purposes.

All **by-chromosome** intermediate results from each step of pipelines are also available on figshare, for celllines 
- **IMR90** in **10kb resolution**
- **GM12878** in **10kb resolution**

## RNAseq processing ([script](rnaseq_processing/pipeline_rnaseq.sh))

### Step 0 Setup
`homedir=/home/nzhou/hic/rao2014/GM12878_10kb/`

### Step 1  Ensembl gene ID and coordinates 

(Optional, only carry out this step if you are using a non-default Ensembl release)

Default release: 90

-  Download Ensembl gene ID and coordinates mapping

To download another release, go to https://www.ensembl.org/biomart/martview/, select "Attributes" on the lefthand column, in the expanded table of "GENE", select "Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)" and "Gene end (bp)"; in the expanded table of "External References", select "NCBI gene ID".
- Save the downloaded file as  [ensembl_map_coding.txt](rnaseq_processing/ensembl_map_coding.txt)


### Step 2  Generate new linear neighborhoods according to Ensembl gene coordinates

Default provided. 

This has nothing to do with cellline, only carry out this step if you are using a non-default Ensembl release

- Designate cutoff

cutoff is the cutoff for linear distance, i.e. if two genes are more than 10000 bps apart then they are not linear neighbors

`cutoff==10000`
- Run script

`python3 genepairs_collapse_linear.py  "$homedir/intra/linear" ./ensembl_map_coding.txt  "$cutoff"`

### Step 3 Write to the genename folder

Default folder is provided at [gene_names](/rnaseq_processing/gene_names/)


### Step 4 Process RNA-seq quantification file

The RNA-seq quantification files are downloaded from ENCODE for each cellline.

Make sure the cell line of your RNAseq file agrees with the cellline for HiC file

[Here](https://www.encodeproject.org/files/ENCFF781YWT/) is an the example RNA-seq quantification file from GM12878.

`python3 prep_rnaseq.py  ./ensembl_map_coding.txt  "$homedir/rnaseq/ENCFF*.tsv"  $homedir/rnaseq/`


### Step 5 Output y as an R matrix

y needs to be in an order that correpsonds to the adjacency matrix from HiC processing

More than one processed RNA-seq file can be supplied at the same time using wildcards

Each RNA-seq file will be read in as a column in a matirx of rnaseq data


```
for chrm in "${chrms[@]}"

do
	Rscript prep_y_for_MRF_intra.R "$chrm" $homedir/intra/  ./gene_names/  $homedir/rnaseq/rnaseq_*.txt
done
```


## HiC processing ([intra script](HiC_processing/pipeline_build.sh))

### Step 0 Setup and download raw data

- Setup home directory
```
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/intra/
cellline=GM12878
resolution=10kb
mkdir $homedir/raw
mkdir $homedir/norm
mkdir $homedir/data
mkdir $homedir/data/y
mkdir $homedir/info
mkdir $homedir/genepairs
mkdir $homedir/linear
```

- Download HiC data from [Rao et al](https://www.ncbi.nlm.nih.gov/pubmed/25497547)

`wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -P $homedir/raw`

You can also download manually from [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525), just save it in `$homedir/raw/`.

- Unzip 

only unzip the desired resolution and filter

`filter=MAPQGE30`

Default read filter is MAPQ>=30. See Supplemental (Extended Experimental Procedure) of Rao et al 2014 (PMID: 25497547) for details

`tar -xvzf $homedir/raw/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -C $homedir/raw/ --strip-components 4 --wildcards GM12878_combined/"$resolution"_resolution_intrachromosomal/chr*/"$filter"/chr*_"$resolution".* `


### Step 1 Normalize HiC contacts and map to genes

use `python3 normalize_and_map_intra.py --help` to check usage and change default arguments such as resolution

`python3 normalize_and_map_intra.py $homedir ../../rnaseq_processing/ensembl_map_coding.txt`


### Step 2 Collapse interaction frequencies by gene pairs using four summary statistics (mean, median, max and min)
Default quantile to cutoff interaction frequency is 0.9. 

i.e. if mean (or median, max, min) interaction frequency of a gene pair is greater than the 90% quantile of all gene pairs in that chromosome

then an edge is inferred between the two genes

`Rscript genepairs_collapse_intra.R $homedir "$quantile" "$resolution"`
To specify a chromosome, use `Rscript genepairs_collapse_intra.R $homedir "$quantile" "$resolution" '19'`, default is all chromosomes.


### Step 3 (optional) Network info
This step reads HiC gene pairs data into edgelists and output graph related information such as edge count and degree distribution. This is NOT required for the PhiMRF analysis

`method="all"`, method is a variable connected to the summary statistics used. In this step, all are considered

`del="FALSE"`, del is a variable connected to whether to delete the linear neighboring gene pairs in the HiC network

```
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quant"_all_TRUE_"$del"
	Rscript get_mat_info_intra.R "$key" $homedir  ../../rnaseq_processing/gene_names/
done
```


### Step 4 Output adjacency matrices from the graphs

Read HiC gene pairs data into edgelists and output adjacency matrix for each gene network

`method="mean"`, (alternatives: `method="median"`, `method="min"`, `method="max"`)

`del="FALSE"`, del=FALSE: Not deleting all linear neighbors from HiC neighbors

```
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quant"_"$method"_TRUE_"$del"
	Rscript prep_mat_for_MRF_intra.R $key  $homedir ../../rnaseq_processing/gene_names 
done
```


## Running PhiMRF
The [vignettes](http://htmlpreview.github.io/?https://github.com/ashleyzhou972/PhiMRF/blob/master/vignettes/Introduction-PhiMRF.html) in the [PhiMRF package](https://github.com/ashleyzhou972/PhiMRF) has very detailed examples of how to use the package. 

Two datasets are required,

- observed count `y`, in our case the processed RNA-seq count, available as `$homedir/intra/data/y/chrm*_y.rds`
- Adjacency matrix of network structure, in our case the processed HiC gene network, available as `$homedir/intra/data/chrm1_0.9_mean_TRUE_FALSE_neighbors_trans.rds`

Example:
`Rscript run/examples.R`


