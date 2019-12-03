###########################################
# For linear neighbors
# read in gene pairs and output adjacency matrix
###########################################
######################################
rm(list=ls())
source('./__get_edge_list.R')
chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')
##update 20190517
##remove '1','7','17','X' from list of chromosomes to loop
##because already done
#new_chrms = chrms[-c(1,7,17,23)]
for (chrm in chrms){
	edgelist = read_edge_list_from_genepairs('./by_gene/cutoff_10kb', chrm)
	mat = output_adjacency_matrix(edgelist, TRUE, chrm)
	savedata('./data/cutoff_10kb/',mat,'../rnaseq_copy.txt', paste0('chrm', chrm))
}
