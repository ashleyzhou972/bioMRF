# Running PhiMRF

All of the following steps are in [`pipeline_PhiMRF.sh`](pipeline_PhiMRF.sh). 

Please read through each step and customize, before running

`bash pipeline_PhiMRF.sh`


The [vignettes](http://htmlpreview.github.io/?https://github.com/ashleyzhou972/PhiMRF/blob/master/vignettes/Introduction-PhiMRF.html) in the [PhiMRF package](https://github.com/ashleyzhou972/PhiMRF) has very detailed examples of how to use the package. 

Two datasets are required,

- Adjacency matrix of network structure: 

    - the processed HiC gene network

    - Example: `$homedir/intra/data/chrm1_0.9_mean_TRUE_FALSE_neighbors_trans.rds`
- Observed count `y`:

    - processed RNA-seq count, 

    - Example: `$homedir/intra/data/y/chrm1_y.rds`

### Designate home dirctory

The same directory as the [set up](../0setup/) step.

Example:
```
homedir=/home/nzhou/hic/rao2014/GM12878_10kb/
```

### Designate quantile and method

```
quantile=0.9
method="mean"
```

### Create results folders

```
mkdir $homedir/intra/results/
mkdir $homedir/inter/results
```

### Before you run PhiMRF in a loop...

1. Please read the [vignette](http://htmlpreview.github.io/?https://github.com/ashleyzhou972/PhiMRF/blob/master/vignettes/Introduction-PhiMRF.html) in the [PhiMRF package](https://github.com/ashleyzhou972/PhiMRF)

2. Please edit the R scripts `run_PhiMRF_*.r` according to documentation in PhiMRF.  

Arguments to edit include: 

- number of total iteration of MCMC

- number of burn-in iterations

- variance of random walk chains

- etc (see documentation for the function `PhiMRF::pmrf()`)

The R scripts to edit:

- [`run_PhiMRF_intra.R`](run_PhiMRF_intra.R)

- [`run_PhiMRF_linear.R`](run_PhiMRF_linear.R)

- [`run_PhiMRF_inter.R`](run_PhiMRF_inter.R)

3. Each run of PhiMRF for each chromosome (pair) could take minutes to hours, depending on the size of the gene network, and the number of iterations. Please tune the model and allow ample time, before committing to a big loop like the following.
 

### Run in loop
```
for chrmA in "${chrms[@]}"
do
	Rscript run_PhiMRF_intra.R $homedir/intra "$chrmA" "$quantile" "$method"
	Rscript run_PhiMRF_linear.R $homedir/intra "$chrmA" 
	for chrmB in "${chrms[@]}"
	do
		Rscript run_PhiMRF_inter.R $homedir/inter "$chrmA" "$chrmB" "$quantile" "$method"
	done
done
```
