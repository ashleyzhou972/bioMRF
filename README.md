# bioMRF: Exploring the Spatial Dependency of Gene Expression
This repository contains 
- scripts to *organize and clean* HiC and RNAseq data
- scripts to run the [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) model on RNAseq and HiC data
- The source code of the R package [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) is available as a separate repo.


## Structure

The pipeline can be separated into three parts

- [set up](0setup/)

- [processing RNAseq data](1rnaseq_processing/)

- [processing HiC data](2hic_processing/)

- [running PhiMRF](3run_PhiMRF/)

The four pipelines are explained in their separate folders, and should be run sequentially.

All intermediate results from each step of pipelines are also available on figshare, for celllines 
- **IMR90** in **10kb resolution**
- **GM12878** in **10kb resolution**


