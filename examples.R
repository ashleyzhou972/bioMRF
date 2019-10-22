####################################################
# An example usage for the PhiMRF model with simulated data
# Naihui Zhou (ashley.n.zhou@gmail.com)
# Last updated 20191021
####################################################
library(PhiMRF)
# irregular_lattice is provided in the package, as an example network structure
adj_trans = transform_small(irregular_lattice)
# Get data size (number of locations)
n = dim(irregular_lattice)[1]
# Simulate y with given parameter values
y = simulate_y(n = n, adj_mat = adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 500)
# Get the parameter space of eta
bounds_e = get_eta_param_space_small(adj_trans)
# Run the pmrf MCMC
ret = pmrf(total_iter = 1000, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
# Make sure the jump counts are within a reasonable range
get_jump_frequency(ret, 1000, 1000)
# Delete initial 200 iterations
new_ret = delete_burn_in(ret, 200)
# Get parameter estimates
print_param_estimates(ret)
