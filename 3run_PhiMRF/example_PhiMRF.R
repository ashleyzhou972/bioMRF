####################################################
# An example usage for the PhiMRF model with Chromosome 1
# Naihui Zhou (ashley.n.zhou@gmail.com)
# Last updated 20191203
####################################################
rm(list=ls())
library(PhiMRF)
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two arguments must be supplied\n", call.=FALSE)
}

root = args[1]
chrm = args[2]
adj_trans = readRDS(paste0(root, '/data/', chrm, '_0.9_mean_TRUE_FALSE_neighbors_trans.rds'))
# Get data size (number of locations)
n = dim(adj_trans)[1]
# Simulate y with given parameter values
y = readRDS(paste0(root, '/data/y/', chrm, '_y.rds'))
# Get the parameter space of eta
bounds_e = get_eta_param_space_small(adj_trans)
# Run the pmrf MCMC
ret = pmrf(total_iter = 20, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
# Make sure the jump counts are within a reasonable range
get_jump_frequency(ret, 20, n)
# Delete initial 200 iterations
new_ret = delete_burn_in(ret, 0)
# Get parameter estimates
print_param_estimates(ret)

saveRDS(ret, "/results/ret_", chrm, ".rds")
