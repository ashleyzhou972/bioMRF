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

## RNAseq processing

### Step 1  Ensembl gene ID and coordinates 

(Optional, only carry out this step if you are using a non-default Ensembl release)

Default release: 90

-  Download Ensembl gene ID and coordinates mapping

To download another release, go to https://www.ensembl.org/biomart/martview/, select "Attributes" on the lefthand column, in the expanded table of "GENE", select "Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)" and "Gene end (bp)"; in the expanded table of "External References", select "NCBI gene ID".
- Save the downloaded file as  [ensembl_map_coding.txt](rnaseq_processing/ensembl_map_coding.txt)


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


### HiC processing
[Rao et al.2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547)

### Running PhiMRF
