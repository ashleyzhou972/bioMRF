# Download raw data and set up home folder

All of the following steps are in [`pipeline_setup.sh`](pipeline_setup.sh). 

Please read through each step and customize, before running 

`bash pipeline_setup.sh`



### Designate the cell line and HiC resolution that you want to run PhiMRF on.

Example: 

```
cellline=GM12878
resolution=10kb
```

### Designate home directory

Example:

`homedir=/home/nzhou/hic/rao2014/GM12878_10kb/`


#### Make folders for RNA-seq data
```
mkdir $homedir/rnaseq
mkdir $homedir/rnaseq/gene_names
```

#### Make folders for intra-chromosomal and inter-chromosomal
```
for int in 'intra', 'inter'
do
	mkdir $homedir/$int
	mkdir $homedir/$int/raw
	mkdir $homedir/$int/norm
	mkdir $homedir/$int/data
	mkdir $homedir/$int/data/y
	mkdir $homedir/$int/info
	mkdir $homedir/$int/genepairs
done
mkdir $homedir/intra/linear
```


### Download intra-chromosomal HiC data (large file)
```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -P $homedir/intra/raw
```
You can also download manually from their [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525). Make sure to save the downloaded file in `$homedir/intra/raw`


### Extract intra-chromosomal HiC data (could take a while)

Only unzip the desired resolution, normalization method and filter.

Example:
```
filter=MAPQGE30
norm=KR
```
See Supplemental 1 (Extended Experimental Procedure) of [Rao et al 2014 (PMID: 25497547)](https://www.ncbi.nlm.nih.gov/pubmed/25497547) for details.

```
tar -xvzf $homedir/intra/raw/GSE63525_"$cellline"_combined_intrachromosomal_contact_matrices.tar.gz -C $homedir/intra/raw/ --strip-components 4 --wildcards GM12878_combined/"$resolution"_resolution_intrachromosomal/chr*/"$filter"/*.{RAWobserved,"$norm"norm}
```

### Download inter-chromosomal HiC data (large file)
```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_"$cellline"_combined_interchromosomal_contact_matrices.tar.gz -P $homedir/inter/raw
```
You can also download manually from their [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525). Make sure to save the downloaded file in `$homedir/inter/raw`

### Extract inter-chromosomal HiC data (could take a while)
```
tar -xvzf $homedir/inter/raw/GSE63525_"$cellline"_combined_interchromosomal_contact_matrices.tar.gz -C $homedir/inter/raw/ --strip-components 4 --wildcards GM12878_combined_interchromosomal/"$resolution"_resolution_interchromosomal/chr*_chr*/"$filter"/*.{RAWobserved,"$norm"norm}
```


### (Optional!) Download Ensembl gene ID and coordinates

Default Ensembl release: 90

This step is optional, only carry out this step if you are want to use another Ensembl release.

To download another release, 

1. Go to [biomart](https://www.ensembl.org/biomart/martview/).
2. Select "Attributes" on the lefthand column, in the expanded table of "GENE", select "Gene stable ID", "Chromosome/scaffold name", "Gene start (bp)" and "Gene end (bp)"; in the expanded table of "External References", select "NCBI gene ID".
3. Save the downloaded file as  [`../1rnaseq_processing/ensembl_map_coding.txt`](../1rnaseq_processing/ensembl_map_coding.txt)


### Download RNA-seq data 

The RNA-seq quantification data are downloaded from ENCODE. Make sure the cell line of your RNA-seq file agrees with the cell line of your HiC data.

If you wish to use non-ENCODE RNA-seq quantification file, make sure the file follows [RSEM's output format](https://www.encodeproject.org/documents/0c78ea4b-9392-421b-a6f3-6c858b6002aa/@@download/attachment/RSEM_Documentation.pdf).

**Note**: The first column should be "gene_id" and second column should be "transcript_id", contrary to the above document.

The quantification metric used in this project is the 8th column: "posterior_mean_count"


Example for GM12878 (two isogenic replicates):
```
wget https://www.encodeproject.org/files/ENCFF781YWT/@@download/ENCFF781YWT.tsv -P $homedir/rnaseq/
wget https://www.encodeproject.org/files/ENCFF680ZFZ/@@download/ENCFF680ZFZ.tsv -P $homedir/rnaseq/
```

### Go to the [next step](../1rnaseq_processing)
