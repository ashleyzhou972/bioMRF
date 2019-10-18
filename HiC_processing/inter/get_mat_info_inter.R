#!/usr/bin/env Rscript
######################################
# For general interHiC networks
########################################
rm(list=ls())
source('../get_edge_list.R')
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  print(args[1])
}
##update 20190517
##remove '1','7','17','X' from list of chromosomes to loop
##because already done
#new_chrms = chrms[-c(1,7,17,23)]

#args should be in the format of "$chrm"_"$quant"_"$method"_TRUE_"$del"

key = args[1]
#key = 'chrm1_chrm7_0.9_all_TRUE_FALSE'
#cat(key, "\n")
keys = strsplit(key, '_')[[1]]
chrm_long1 = keys[1] #chrm_long contains the string "chrm"
chrm_long2 = keys[2]
chrm1 = substr(chrm_long1, 5,6) #chrm only contains the number of the chromosome
chrm2 = substr(chrm_long2, 5,6)
chrms = c(chrm1, chrm2)
quant = as.numeric(keys[3])
method = keys[4]
add_singles = as.logical(keys[5])
del_neighbors = as.logical(keys[6])

homefolder = '/home/nzhou/hic/rao2014/'
output_folder = './IMR90_10kb/inter/by_gene/info/'
genename_folder = './gene_names/'
if (compare_chrms(chrm1, chrm2)){

	info.ecount = c()
	info.degrees_median = c()
	info.degrees_mean = c()
	info.degrees_sd = c()
	#info.percent_loop = c()

	if (method=='all'){
	  #read in unique loop network
	  #loop_neighbor_edge_list = read.table(paste0('./IMR90_10kb/hiccups/thres10/edge_list_',chrm_long,'.txt'))
	  #loop_net0 = graph_from_edgelist(as.matrix(loop_neighbor_edge_list), directed=F)
	  #loop_net = simplify(loop_net0)

	  edgelist_list = read_edge_list_from_genepairs_inter('./IMR90_10kb/inter/by_gene/genepairs/', chrm1, chrm2, quant, method)
	  method_names = c("mean", "median", "max", "min")
	  i = 1
	  for (el in edgelist_list$el){
	    net = output_graph(homefolder, genename_folder, el, add_singles, del_neighbors, chrms)
	    info.ecount[i] = ecount(net)
	    degrees = degree(net)
	    info.degrees_median[i] = median(degrees)
	    info.degrees_mean[i] = mean(degrees)
	    info.degrees_sd[i] = sd(degrees)
	    #inter_with_loop = intersection(loop_net, net, keep.all.vertices = F)
	    #inter_with_linear = intersection(linear_net, net, keep.all.vertices = F)
	    #cat("\n ratio of edges from loop calling in edges from HiC\n")
	    #r1 = ecount(inter_with_loop)/ ecount(net)
	    #cat(r1)
	    #cat("\n ratio of edges from HiC in edges from loop calling\n")
	    #r2 = ecount(inter_with_loop)/ecount(loop_net)
	    #cat(r2)
	    #cat("\n ratio of edges from linear neighbors in edges from HiC\n")
	    #r3 = ecount(inter_with_linear)/ ecount(net)
	    #cat(r3)
	    #cat("\n ratio of edges from HiC in edges from linear neighbors\n")
	    #r4 = ecount(inter_with_linear)/ ecount(linear_net)
	    #cat(r4)
	    #cat("\n")
	    #info.percent_loop[i] = r2
	    #info.percent_linear[i] = r4
	    i = i + 1
	  }  
	} else {
	  stop("Use 'all' as method in input key\n")
	}

	#table = cbind(info.ecount, info.degrees_median, info.degrees_mean, info.degrees_sd, 
	#              info.percent_loop, info.percent_linear)
	table = cbind(info.ecount, info.degrees_median, info.degrees_mean, info.degrees_sd)
	name = paste0('chrm', chrm1, '_', 'chrm', chrm2, '_', quant, '_', method,'_', add_singles, '_', del_neighbors)
	saveRDS(table, paste0(output_folder, '/', name, '_info.rds'))
} else {
	cat("Invalid chromosome pair,", chrm1, chrm2, "\n")
}
