rm(list=ls())
library(Matrix, warn.conflict=F)
library(PhiMRF)
args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("Four arguments must be supplied\n", call.=FALSE)
}

root = args[1]
chrm = args[2]
quant = args[3]
method = args[4]


#########################################################
# CHANGE THESE!
Ti = 5000 # total iteration
vars = c(1, 0.5, NA, 2) # variance for w, alpha, eta and tau2
B = 1000 # burn in
########################################################


adj_trans = readRDS(paste0(root, '/data/', chrm, '_',quant, '_', method, '_TRUE_FALSE_neighbors_trans.rds'))
# Get data size (number of locations)
n = dim(adj_trans)[1]
# Simulate y with given parameter values
y = readRDS(paste0(root, '/data/y/', chrm, '_y.rds'))
# Get the parameter space of eta
bounds_e = get_eta_param_space_small(adj_trans)
# Run the pmrf MCMC
ptm <- proc.time()
ret = pmrf(total_iter = Ti, N = n , y, adj_trans, vars, bounds_e)
cat("Time used:\n")
print(proc.time() - ptm)
# Make sure the jump counts are within a reasonable range
get_jump_frequency(ret, Ti, n)
# Delete initial 200 iterations
new_ret = delete_burn_in(ret, B)
# Get parameter estimates
cat("Parameter Estimate results:\n")
print_param_estimates(ret, B)

saveRDS(ret, paste0(root, "/results/ret_intra_", chrm, ".rds"))
