rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Three arguments must be supplied\n", call.=FALSE)
} 
source('./__subset_mat.R')

chrms = as.character(seq(1,22,1))
chrms = c(chrms, 'X')

TADfolder = args[1]
by_chrm_data_folder = args[2] # data path for full chromosomes gene networks
outfolder = args[3]


#TADs: intra-only

get_dimname_1<-function(matrix){
  return(attributes(matrix)$Dimnames[[1]])
}

for (chrm in chrms){
  TAD_mats = list()
  TAD_y = c()
  chrm_long = paste0("chrm", chrm)
  key = paste0(chrm_long, '_0.9_mean_TRUE_FALSE')
  TADgenes = read.table(paste0(TADfolder, '/', chrm_long, '_TAD.txt'), header = F)
  # get subset of each TAD and then concatenate together
  # as block diagonal matrix
  # then get index, output y 
  colnames(TADgenes)<-c("gene", "domain")
  gene_domains= aggregate(gene~domain, TADgenes, paste, collapse = "--")
  genes = as.character(TADgenes[,1])
  #TADsub = induced_subgraph(net, vids=genes)
  #group1 = gene_domains[1,2]
  #vec1 = unlist(strsplit(group1, split="--", fixed=T))
  #TADsub_intra = induced_subgraph(TADsub, vec1)
  for (group in gene_domains[,2]){ 
    vec = unlist(strsplit(group, split="--", fixed=T))
    #print(vec)
    ret = subset_mat(by_chrm_data_folder, key, chrm_long, TRUE, vec, NULL, NULL, NULL )
    #sub_net = induced_subgraph(TADsub, vec)
    #TADsub_intra = union(TADsub_intra, sub_net)
    TAD_mats = c(TAD_mats, ret$net)
    TAD_y = rbind(TAD_y, ret$y)
  }
  big = bdiag(TAD_mats) #The small matrices are already transformed
  names = unlist(sapply(TAD_mats, get_dimname_1, simplify =T))
  attributes(big)$Dimnames[[1]]<-names
  attributes(big)$Dimnames[[2]]<-names
  #But since they are block diagonal, the transformation still holds.
  #check if dimensions agree
  if (dim(big)[1]!=length(genes) || dim(big)[2]!=length(genes)){
    stop("Dimension of intra-TAD matrix does not agree with number of TAD genes\n")
  }
  #check if y and mat have the same order
  if (sum(attributes(big)$Dimnames[[1]]!=rownames(TAD_y))>0){
    stop("Gene names of resulting matrix not the same as y\n")
  }
  #save
  saveRDS(big, file = paste0(outfolder, "/", paste0("TADs_intra_", chrm_long), "_neighbors_trans.rds"))
  saveRDS(TAD_y, file = paste0(outfolder, "/y/", paste0("TADs_intra_", chrm_long), "_y.rds"))
}
