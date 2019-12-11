# Process HiC data

All of the following steps are in [`pipeline_hic_intra.sh`](intra/pipeline_hic_intra.sh) or [`pipeline_hic_inter.sh`](inter/pipeline_hic_inter.sh).

Please read through each step and customize, before running

```
bash pipeline_hic_intra.sh`
```

or 

```
bash pipeline_hic_intra.sh`
```

The following uses **intra** as example.


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

### Designate resolution

Example:
```
resolution=10kb
```


## Normalize HiC contacts and map to genes

Use `python3 normalize_and_map_intra.py --help` to check usage.

Example:
```
python3 normalize_and_map_intra.py $homedir/intra $rnaseq --resolution "$resolution" --norm "KR" --overlap 0.1
```


## Collapse interaction frequencies by gene pairs using four summary statistics (mean, median, max and min)

1. Set quantile

```
quantile=0.9
```
Default quantile to cutoff interaction frequency is 0.9. 

i.e. if mean (or median, max, min) interaction frequency of a gene pair is greater than the 90% quantile of all gene pairs in that chromosome, then an edge is inferred between the two genes. 

2. Collapse interactions (could take a while)

```
Rscript genepairs_collapse_intra.R $homedir/intra "$quantile" "$resolution"
```

The above default command runs through all chromosomes. A fourth argument can be supplied for a single chromosome.

Example for the 19th chromosome: `Rscript genepairs_collapse_intra.R $homedir/intra "$quantile" "$resolution" '19'`.



### Generate linear neighborhoods 

This step outputs the adjacency matrices for gene networks based on linear distance on each chromosome, as opposed to 3D HiC distance.

1. Set cutoff

Example:
```
cutoff=10000 
```
Cutoff is the cutoff for linear neighbors, for example, if two genes are more than 10000 bps apart then they are not linear neighbors.

2. Get linear gene pairs 
```
python3 genepairs_collapse_linear.py  "$homedir/intra/linear" $rnaseq  "$cutoff"
```

3. Output linear adjacency matrix
```
Rscript prep_mat_for_MRF_linear.R $homedir/intra/linear $rnaseq $homedir/rnaseq/gene_names/ "TRUE"
```


### (Optional!) Output Network info

This step reads HiC gene pairs data into edgelists and output graph related information such as edge count and degree distribution. This step is not required for the PhiMRF analysis.

1. Set method (summary statistic)

Only "all" is accepted in this step
```
method="all"
```

2. Set del

del is a variable connected to whether to delete the linear neighboring gene pairs in the HiC network. The option to delete linear neighbors in the HiC network isolates the **long-range** gene-gene interactions in the genome.

Example:
```
del="FALSE"
```

3. Run through every chromosome
```
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quantile"_all_TRUE_"$del"
	Rscript get_mat_info_intra.R "$key" $homedir/intra  $homedir/rnaseq/gene_names/
done
```


### Output adjacency matrices

1. Set method (summary statistic)

```
method="mean"
```
alternatives:
```
method="median"
method="min"
method="max"
```

2. Set del

del is a variable connected to whether to delete the linear neighboring gene pairs in the HiC network. The option to delete linear neighbors in the HiC network isolates the **long-range** gene-gene interactions in the genome.

Example:
```
del="FALSE"
```

3. Output adjacency matrix
```
for chrm in "${chrms[@]}"
do
	key="$chrm"_"$quantile"_"$method"_TRUE_"$del"
	Rscript prep_mat_for_MRF_intra.R $key  $homedir/intra $homedir/rnaseq/gene_names 
done
```
