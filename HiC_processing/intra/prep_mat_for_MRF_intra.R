#!/usr/bin/env Rscript
######################################
# For general intra HiC networks
########################################
rm(list=ls())
source('../get_edge_list.R')
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  print(args[1])
}
homefolder = '/home/nzhou/hic/rao2014/'
genename_folder = '/home/nzhou/hic/rao2014/gene_names/'
##update 20190517
##remove '1','7','17','X' from list of chromosomes to loop
##because already done
#new_chrms = chrms[-c(1,7,17,23)]

#args should be in the format of "$chrm"_"$quant"_"$method"_TRUE_"$del"

key = args[1]
#key = 'chrm22_0.9_mean_TRUE_TRUE'
cat(key, "\n")
keys = strsplit(key, '_')[[1]]
chrm = substr(keys[1], 5,6)
quant = as.numeric(keys[2])
method = keys[3]
add_singles = as.logical(keys[4])
del_neighbors = as.logical(keys[5])

output_folder = './IMR90_10kb/intra/by_gene/data/'
  
edgelist_list = read_edge_list_from_genepairs('./IMR90_10kb/intra/by_gene/genepairs/', chrm, quant, method)
if (method=='all'){
  method_names = c("mean", "median", "max", "min")
  i = 1
  for (el in edgelist_list$el){
    net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrm)
    mat = output_adjacency_matrix(genename_folder, net, chrm, add_singles)
    name = paste0('chrm', chrm,  '_', quant, '_', method_names[i],'_', add_singles, '_', del_neighbors)
    #print(name)
    save_mat(output_folder,mat, name)
    i = i + 1
  }  
} else {
  el = edgelist_list$el
  net = output_graph(homefolder,genename_folder, el, add_singles, del_neighbors, chrm)
  mat = output_adjacency_matrix(genename_folder, net, chrm, add_singles)
  name = paste0('chrm', chrm,  '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
  save_mat(output_folder, mat, name)
}

