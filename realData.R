## Paper: Universal inference for model selection on networks
## Author: Eric Yanchenko
## Akita International University

source("~/Documents/Research/network_model_selection/functions.R")

load("~/Documents/Research/Srijan/CP_clean/Data/polblogs_adj.RData")

nrow(A)
sum(A)/(nrow(A)*(nrow(A)-1))

gamma = 0.5

eval_mc(A, "ER", "CL",     2, 2, gamma, nreps = 1000, detectCores()-1)
eval_mc(A, "ER", "SBM",    2, 2, gamma, nreps = 1000, detectCores()-1)
eval_mc(A, "CL", "DCBM",   2, 2, gamma, nreps = 1000, detectCores()-1)
eval_mc(A, "SBM", "DCBM",  2, 2, gamma, nreps = 1000, detectCores()-1)
eval_mc(A, "DCBM", "DCBM", 2, 3, gamma, nreps = 1000, detectCores()-1)
eval_mc(A, "DCBM", "PABM", 2, 2, gamma, nreps = 1000, detectCores()-1)