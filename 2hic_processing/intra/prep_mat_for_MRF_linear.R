###########################################
# For linear neighbors
# read in gene pairs and output adjacency matrix
###########################################
######################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
	  stop("Four arguments must be supplied\n", call.=FALSE)
} 
source('./__get_edge_list.R')
chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')

homedir = args[1]
rnaseq_file = args[2]
genename_folder = args[3]
singleton=as.logical(args[4])

read_edge_list_from_genepairs_linear<-function(infolder, chrm){
  result = as.data.table(read.table(paste0(infolder, '/linear_chr', chrm,'.genepairs'), header = F, sep='\t'))
  edgelist = as.matrix(result)
  return(edgelist)
}

for (chrm in chrms){
  el = read_edge_list_from_genepairs_linear(homedir, chrm)
  net = output_graph(homedir, genename_folder, el, singleton, FALSE, chrm)
  mat = output_adjacency_matrix(genename_folder, net, chrm, singleton)
  name = paste0('chrm', chrm)
  save_mat(homedir, mat, name)

}
