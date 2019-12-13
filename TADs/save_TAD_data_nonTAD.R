rm(list=ls())
source('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/change_neighbor_mat.R')
datafolder = '/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/data/'
#chrms = as.character(seq(1,22,1))
#chrms = c(chrms, 'X')
chrms = c("1")

#nonTADs
TADfolder = '/home/nzhou/hic/rao2014/IMR90_10kb/hiccups/TADs/TADgenes/'
for (chrm in chrms){
  chrm_long = paste0("chrm", chrm)
  nonTADgenes = read.table(paste0(TADfolder, '/', chrm_long, '_nonTAD.txt'), header = F)[,1]
  ret = subset_mat(datafolder, paste0(chrm_long, '_0.9_mean_TRUE_FALSE'), chrm_long, TRUE, nonTADgenes, 
                   "/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/nonTAD_data", 
                   paste0("nonTADs_", chrm_long) , paste0("nonTADs_", chrm_long))
  
}