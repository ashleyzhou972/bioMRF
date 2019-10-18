key='chrm22'
chrm_long = key
#cat(key, "\n")
################################################
# Run nonTAD  genes in each chromosome and compare
################################################
source('/home/nzhou/hic/IMR90/work/MRF_HIC_GE/analysis/post_analysis.R')
#install.packages("~/hic/R_dev/bioMRF_0.1.0.tar.gz", repos = NULL, type = "source")
library(PhiMRF)
set.seed(234)

datafolder = "/home/nzhou/hic/rao2014/IMR90_10kb/intra/by_gene/nonTAD_data/"
new_neighbor = readRDS(paste0(datafolder, "/", "nonTADs_", chrm_long, "_neighbors_trans.rds"))
y = readRDS(paste0(datafolder,"/y/nonTADs_", chrm_long, "_y.rds"))
n = nrow(y)

nnz = length(attributes(new_neighbor)$x)
cat("Number of nonzero elements in neighborhood matrix is ", nnz, "\n")

#bounds_e = PhiMRF::get_eta_param_space_small(new_neighbor)

new_new_neighbor = preprocess_big(new_neighbor)$new_m
#########toggle#####################
vars<-c(2, 1.5, 0.02, 2)
total_iter = 2000
B = 40
shouldsave = F
#savepath = paste0('./returns/nonTADs/ret_',key,'.rds')
##########end toggle###############

#there are four steps of metropolis in each iteration (including double metropolis)

#parameters for the prior distributions (uniform), for alpha, eta and tau


#inital guess for alpha, beta and tau^2
inis1 = c(0.1,0.0,0.1)
#saverds(inis1, file=paste('../results/',date,'/tests/','/inis_os1',counter,'.rds', sep = ''))
winis = rnorm(n, inis1[1], sqrt(inis1[3]))
ret <- pmrf(total_iter, n, y, new_neighbor, vars, bounds_e, vars = vars)
ret1 <- pmrf(total_iter, n, y, new_neighbor, vars, bounds_e)

cat('###########################################\n')
cat("random walk variance:", vars, "\n")
cat('###########################################\n')
jump_count = get_jump_frequency(ret, total_iter, n)
print(jump_count)

cat('###########################################\n')
cat('Burn-in:', B, '\n')

bioMRF::print_param_estimates(ret, burn_in = B)


if (shouldsave) {
	  saveRDS(ret,savepath)
}
