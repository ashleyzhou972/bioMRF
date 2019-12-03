# bioMRF: Exploring the Spatial Dependency of Gene Expression

- This repository contains scripts to *organize and clean* HiC interaction data.
- The R package [PhiMRF](https://github.com/ashleyzhou972/PhiMRF) is used to do the actual calculation of the SIE. 
- This repository contains examples and demos for using the [PhiMRF package](https://github.com/ashleyzhou972/PhiMRF) on RNA-seq and HiC data.

Each steps below are done separately for intra-chromosomal cases and inter-chromosomal cases.

### Map HiC raw contact matrices into spatial gene networks (intra and inter)

- Download HiC contact matrix data.

In this study, we use raw HiC contact matrices from [Rao *et al.*2014](https://www.ncbi.nlm.nih.gov/pubmed/25497547), in which HiC data of a variety of celllines and resolutions are available. We use the IMR90 cell line and 10kb resolution as example. 

The choice of cell line depends on the availability of RNA-seq data of the same cell line. We recommend higher resolutions to ensure the fine mapping of genes to HiC contact bins.

### 
