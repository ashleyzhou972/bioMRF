rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two argument must be supplied\n", call.=FALSE)
} else if (length(args)==2) {
  key = args[1]
  homefolder = args[2]
}
#key = "GO:0007165_False"
#anno = read.table(paste0('./', key,'.tsv'), sep = "\t", header = T)
#signal_genes = unique(anno$GENE.PRODUCT.ID)
#write.table(signal_genes, paste0("./",key,"_genes_uniprot.txt"), row.names = F, col.names = F, quote=F)
#### use online mapping tool on uniprot.org

source('../TADs/__subset_mat.R')
datafolder = paste0(homefolder, "/functional/data/")
signalFolder = paste0(homefolder, "functional/annotations/ensembl_list/")
##read signal genes

signal_genes = as.character(read.table(paste0(signalFolder, '/', key,'_genes_ensembl.txt'), header = T)[,1] )

key_all = "ALL"
savefolder = datafolder
savepath = paste0(savefolder, key, "_data/")
if (!dir.exists(savepath)){
  dir.create(savepath, showWarnings = T)
  dir.create(paste0(savepath, '/y/'), showWarnings=T)
  cat("directory created\n")
  ret1 = subset_mat(datafolder, paste0(key_all, '_0.9_mean_TRUE_FALSE'), key_all, FALSE, signal_genes, 
                    savepath, key, key)
} else {
  #cat("neighbors_trans file already exists\n")
  ret1 = subset_mat(datafolder, paste0(key_all, '_0.9_mean_TRUE_FALSE'), key_all, FALSE, signal_genes, 
                    NULL, NULL, NULL)
}

#save connectivity stats: vcount,ecount, ecount/max_e, ecount/vcount, 
#median_degree, median_degree/vcount, num_singleton
#max_degree, max_degree/vcount
library(igraph, warn.conflicts = F)
graph = graph_from_adjacency_matrix(ret1$net, mode = "undirected", weighted = T)

vcount = vcount(graph)
max_e = choose(vcount,2)
ecount = ecount(graph)
e_ratio = round(ecount/max_e,3)
ev_ratio = round(ecount/vcount,3)
med_degree = median(degree(graph))
max_degree = max(degree(graph))
med_deg_ratio = round(med_degree/vcount,3)
max_deg_ratio = round(max_degree/vcount,3)
num_sing = sum(degree(graph)==0)

cat(paste(vcount, ecount, e_ratio, ev_ratio, med_degree, med_deg_ratio, max_degree, max_deg_ratio, num_sing, "\n", sep = "\t"))
