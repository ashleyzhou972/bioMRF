rm(list=ls())
source('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/change_neighbor_mat.R')
datafolder = '/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/data/'
#chrms = as.character(seq(1,22,1))
#chrms = c(chrms, 'X')
chrms = c("1")

#TADs
TADfolder = '/home/nzhou/hic/rao2014/IMR90_10kb/hiccups/TADs/TADgenes/'
for (chrm in chrms){
  chrm_long = paste0("chrm", chrm)
  key = paste0(chrm_long, '_0.9_mean_TRUE_FALSE')
  TADgenes = read.table(paste0(TADfolder, '/', chrm_long, '_TAD.txt'), header = F)[,1]
  ret = subset_mat(datafolder, key, chrm_long, TRUE, TADgenes, 
                   "/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/TAD_data/temp/", 
                   paste0("TADs_", chrm_long) , paste0("TADs_", chrm_long))
  
}
