## Running PhiMRF
The [vignettes](http://htmlpreview.github.io/?https://github.com/ashleyzhou972/PhiMRF/blob/master/vignettes/Introduction-PhiMRF.html) in the [PhiMRF package](https://github.com/ashleyzhou972/PhiMRF) has very detailed examples of how to use the package. 

Two datasets are required,

- observed count `y`, in our case the processed RNA-seq count, available as `$homedir/intra/data/y/chrm*_y.rds`
- Adjacency matrix of network structure, in our case the processed HiC gene network, available as `$homedir/intra/data/chrm1_0.9_mean_TRUE_FALSE_neighbors_trans.rds`

Example:
`Rscript run/examples.R`
