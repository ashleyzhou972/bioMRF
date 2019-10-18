####
# Output all gene names
####
rm(list=ls())
setwd("/home/nzhou/hic/rao2014/gene_names/")

chrms = as.character(seq(1,22))
chrms = c(chrms, 'X')

all_names = c()
for (chr in chrms){
  names = readRDS(paste0("./chrm", chr, "_genes.rds"))
  all_names = c(all_names, names)
  #print(length(all_names))
}

write.table(all_names, "./all_gene_names.txt", row.names = F, col.names = F, quote=F, sep='\n')
