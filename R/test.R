key = 'chrm22'
chrm_long = key
cat(key, "\n")
################################################
# Run TAD  genes in each chromosome and compare
################################################
#library(Matrix)
#setwd('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/real_by_chrm')
#source('./simulation_1000.R')
source('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/post_analysis.R')
set.seed(234)

datafolder = "/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/TAD_intra_data/"
new_neighbor = readRDS(paste0(datafolder, "/", "TADs_intra_", chrm_long, "_neighbors_trans.rds"))
y = readRDS(paste0(datafolder,"/y/TADs_intra_", chrm_long, "_y.rds"))
n = nrow(y)

nnz = length(attributes(new_neighbor)$x)
cat("Number of nonzero elements in neighborhood matrix is ", nnz, "\n")

evalues = eigen(new_neighbor)$values #do eigen using the sparse matrix format
eta_l = 1/min(evalues)
cat('lower bound for eta', eta_l, '\n')
eta_u = 1/max(evalues)
cat('upper bound for eta', eta_u, '\n')
load_y = TRUE

#########toggle#####################
vars<-c(2, 1.5, 0.02, 2)
total_iter = 2000
B = 400
shouldsave = F
#savepath = paste0('./returns/TADs_intra/ret_',key,'.rds')
##########end toggle###############

#there are four steps of metropolis in each iteration (including double metropolis)

#parameters for the prior distributions (uniform), for alpha, eta and tau
bounds_a = c(-10,10)
bounds_e = c(eta_l, eta_u)
bounds_t = c(0,10) #This is the range for tau

#inital guess for alpha, beta and tau^2
inis1 = c(0.1,0.0,0.1)
#saverds(inis1, file=paste('../results/',date,'/tests/','/inis_os1',counter,'.rds', sep = ''))
winis = rnorm(n, inis1[1], sqrt(inis1[3]))
ret <- bioMRF::dm_call_wrapper(total_iter, n, y, new_neighbor, vars, bounds_a, bounds_e, bounds_t, inis1, winis)

cat('###########################################\n')
cat("random walk variance:", vars, "\n")
cat('###########################################\n')
jump_count = get_jump_frequency(ret, total_iter, n)
print(jump_count)

cat('###########################################\n')
cat('Burn-in:', B, '\n')
newret <-delete_burn_in(ret, B, new_neighbor)

print_param_estimates(newret)


if (shouldsave) {
	  saveRDS(ret,savepath)
}
