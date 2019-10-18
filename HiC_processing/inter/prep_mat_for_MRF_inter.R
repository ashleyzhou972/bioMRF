#!/usr/bin/env Rscript
######################################
# For general inter HiC networks
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
chrm1 = substr(keys[1], 5,6)
chrm2 = substr(keys[2], 5,6)
quant = as.numeric(keys[3])
method = keys[4]
add_singles = as.logical(keys[5])
del_neighbors = as.logical(keys[6])

output_folder = './IMR90_10kb/inter/by_gene/data/'
if ( compare_chrms(chrm1, chrm2)) {
	edgelist_list = read_edge_list_from_genepairs_inter('./IMR90_10kb/inter/by_gene/genepairs/', chrm1, chrm2, quant, method)
	chrms = c(chrm1, chrm2)
	if (method=='all'){
	  method_names = c("mean", "median", "max", "min")
	  i = 1
	  for (el in edgelist_list$el){
	    net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrms)
	    mat = output_adjacency_matrix(genename_folder, net, chrms, add_singles)
	    name = paste0('chrm', chrm1, '_', 'chrm', chrm2, '_', quant, '_', method_names[i],'_', add_singles, '_', del_neighbors)
	    #print(name)
	    save_mat(output_folder,mat, name)
	    i = i + 1
	  }  
	} else {
	  el = edgelist_list$el
	  net = output_graph(homefolder,genename_folder, el, add_singles, del_neighbors, chrms)
	  mat = output_adjacency_matrix(genename_folder, net, chrms, add_singles)
	  name = paste0('chrm', chrm1, '_', 'chrm', chrm2, '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
	  save_mat(output_folder, mat, name)
	}
} else {
	cat("Invalid chromosome pair,", chrm1, chrm2, "\n")
}
