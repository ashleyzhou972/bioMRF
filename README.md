# bioMRF: Spatial Dependency of Quantitative Gene Expression in the 3D Genome

### This repository contains 
- scripts to organize and clean HiC and RNAseq data
- scripts to run [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) on RNAseq and HiC data

The R package [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) is the statistical model we developed to detect spatial dependency in count data observed on spatial structures. The source code of the package as well as documentation is available as a separate [repository](https://github.com/ashleyzhou972/PhiMRF).

### Prerequisites: 

- [PhiMRF](https://github.com/ashleyzhou972/PhiMRF)
- python3
- R
- bash

### This pipeline has four parts:

- [Set up](0setup/)

- [Processing RNA-seq data](1rnaseq_processing/)

- [Processing HiC data](2hic_processing/)

- [Running PhiMRF](3run_PhiMRF/)

The four parts are explained in their separate folders, and should be run sequentially.

All intermediate results from each step of this pipeline are available on figshare, for the following celllines.

- **IMR90** in **10kb resolution**

- **GM12878** in **10kb resolution**

### References

- RNAseq data comes from [ENCODE](https://www.encodeproject.org)

ENCODE Project Consortium. "An integrated encyclopedia of DNA elements in the human genome." Nature 489.7414 (2012): 57.

- HiC data and normalization method comes from [Rao et al, 2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547)

Rao, Suhas SP, et al. "A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping." Cell 159.7 (2014): 1665-1680.
