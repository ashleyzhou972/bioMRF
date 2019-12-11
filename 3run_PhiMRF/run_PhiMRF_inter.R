rm(list=ls())
source('../2hic_processing/inter/__get_edge_list.R')
library(Matrix, warn.conflict=F)
library(PhiMRF)
args = commandArgs(trailingOnly=TRUE)

if (length(args)<4) {
  stop("Four arguments must be supplied\n", call.=FALSE)
}

root = args[1]
chrmA = args[2]
chrmB = args[3]
quant = args[4]
method = args[5]

#########################################################
# CHANGE THESE!
Ti = 20 # total iteration
vars = c(1, 0.5, NA, 2) # variance for w, alpha, eta and tau2
B = 5 # burn in
########################################################

if (compare_chrms(chrmA, chrmB)){
  adj_trans = readRDS(paste0(root, '/data/', chrmA, '_', chrmB, '_',quant, '_', method, '_TRUE_FALSE_neighbors_trans.rds'))
  # Get data size (number of locations)
  n = dim(adj_trans)[1]
  # Simulate y with given parameter values
  y = readRDS(paste0(root, '/data/y/', chrmA, '_', chrmB, '_y.rds'))
  # Get the parameter space of eta
  bounds_e = get_eta_param_space_small(adj_trans)
  # Run the pmrf MCMC
  ptm <- proc.time()
  ret = pmrf(total_iter = Ti, N = n , y, adj_trans, vars, bounds_e)
  cat("Time used:\n")
  print(proc.time() - ptm)
  # Make sure the jump counts are within a reasonable range
  cat("Jump frequencies:\n")
  print(get_jump_frequency(ret, Ti, n))
  # Delete initial 200 iterations
  new_ret = delete_burn_in(ret, B)
  # Get parameter estimates
  cat("Parameter Estimate results:\n")
  print_param_estimates(ret, B)
  
  saveRDS(ret, paste0(root, "/results/ret_inter_", chrmA, '_', chrmB , ".rds"))
}