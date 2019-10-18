rm(list=ls())
TADfolder = '/home/nzhou/hic/rao2014/IMR90_10kb/hiccups/TADs/TADgenes/'
chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')
#chrm = '1'
for (chrm in chrms){
  chrm_long = paste0("chrm", chrm)
  TADgenes = read.table(paste0(TADfolder, '/', chrm_long, '_TAD.txt'), header = F)
  # get subset of each TAD and then concatenate together
  # as block diagonal matrix
  # then get index, output y 
  colnames(TADgenes)<-c("gene", "domain")
  gene_domains= aggregate(gene~domain, TADgenes, paste, collapse = "--")
}