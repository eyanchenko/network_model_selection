## Data fission and UI for model selection on networks
## Eric Yanchenko
## May 8, 2026

#' @title Generate networks from models
#' @description Generates network from various random-graph models
#' @param model random-graph model
#' @param params (list of) parameters corresponding to the particular model
#' @return adjacency matrix
#' @export
generateA <- function(model=c("ER", "CL", "SBM", "DCBM", "PABM"), params){
  
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name: Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (=DCBM),
    popularity adjusted block model (=PABM)")
  }
  
  if(model=="ER"){
    n = params$n
    p = params$p
    P = matrix(p, nrow=n, ncol=n)
  }else if (model=="CL"){
    psi = params$psi
    n = length(psi)
    P = as.vector(psi %*% t(psi))
  }else if(model=="SBM"){
    B = params$B
    C = params$C
    n = length(C)
    
    K = length(unique(C))
    C <- as.factor(C)
    
    # Membership matrix
    Z <- matrix(model.matrix(~C-1), nrow=n, ncol=K)
    
    P <- Z%*%B%*%t(Z)
  }else if(model=="DCBM"){
    psi = params$psi
    B   = params$B
    C   = params$C
    
    n = length(psi)
    K = length(unique(C))
    C <- as.factor(C)
    
    # Membership matrix
    Z <- matrix(model.matrix(~C-1), nrow=n, ncol=K)
    
    P <- Z%*%B%*%t(Z)
    P <- psi%*%t(psi) * P
  }else if(model=="PABM"){
    lambda = params$lambda
    C      = params$C
    n      = nrow(lambda)
    K      = ncol(lambda)
    
    P <- matrix(0, n, n)
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        P[i,j] <- lambda[i, C[j]] * lambda[j, C[i]]
        P[j,i] <- P[i,j]
      }
    }
    
    P = P/max(P)
  }
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  A[lower.tri(A)] <- 0
  A = A + t(A)
  diag(A) <- 0
  return(A)
}

#' @title Log-likelihood
#' @description Computes the log-likelihood for an observed matrix
#' and given random-graph model (and parameters) using the Poisson approximation.
#' @param A adjacency matrix
#' @param model random-graph model
#' @param params (list of) parameters corresponding to the particular model
#' @return evaluation of log-likelihood
#' @export
llike <- function(A, model=c("ER", "CL", "SBM", "DCBM", "PABM"), params){
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name: Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (=DCBM),
    popularity adjusted block model (=PABM)")
  }
  
  
  if(model=="ER"){
    n = params$n
    p = params$p
    P = matrix(p, nrow=n, ncol=n)
  }else if (model=="CL"){
    psi = params$psi
    n = length(psi)
    P = psi %*% t(psi)
  }else if(model=="SBM"){
    B = params$B
    C = params$C
    n = length(C)
    
    K = length(unique(C))
    C <- as.factor(C)
    
    # Membership matrix
    Z <- matrix(model.matrix(~C-1), nrow=n, ncol=K)
    
    P <- Z%*%B%*%t(Z)
  }else if(model=="DCBM"){
    psi = params$psi
    B   = params$B
    C   = params$C
    
    n = length(psi)
    K = length(unique(C))
    C <- as.factor(C)
    
    # Membership matrix
    Z <- matrix(model.matrix(~C-1), nrow=n, ncol=K)
    
    P <- Z%*%B%*%t(Z)
    P <- psi%*%t(psi) * P
  }else if(model=="PABM"){
    C      = params$C
    lambda = params$lambda
    n      = length(C)
    
    P = matrix(0, nrow = n, ncol = n)
    for (i in 1:n){  	# populate diagonal entries
      r = C[i]
      P[i,i] = (1/2)*(lambda[i,r]^2)
      }	# divide by 2 --- see notes section 3.3
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        p = lambda[i,C[j]]*lambda[j,C[i]]
        P[i,j] = p
        P[j,i] = P[i,j]
      }}
     
  }
  
  eps = 1e-6
  P[P < eps] <- eps
  P[P > 1-eps] <- 1-eps
  
  return(
    -sum(P[upper.tri(P)]) + sum(A[upper.tri(A)] * log(P[upper.tri(P)]))
  )
}

# Helper function for PABM estimation. From Sengupta and Chen (2018)
f.PA<-function(A,b){	
  K<-max(b)       # no. of communities
  N<-nrow(A)      # no. of nodes
  M<-matrix(NA,nrow=N,ncol=K)  # popularity matrix
  O<-matrix(NA,nrow=K,ncol=K)  # community interaction matrix
  for (i in 1:N){		# calculate M
    for (r in 1:K){
      nodes = which(b == r)
      M[i,r] = sum(A[i,nodes])
    }}
  for (r in 1:K){		# calculate O
    for (s in r:K){
      nodes1 = which(b == r)
      nodes2 = which(b == s)
      O[r,s] = sum(A[nodes1,nodes2])
      O[s,r] = O[r,s]
    }}
  list(M=M, O=O)}

#' @title Estimate model parameters
#' @description Estimates model parameters for random graph model
#' @param A adjacency matrix
#' @param model random-graph model
#' @param C community labels (if applicable)
#' @return estimated parameter values
#' @export
estparam <- function(A, model=c("ER", "CL", "SBM", "DCBM", "PABM"), C){
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name: Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (=DCBM),
    popularity adjusted block model (=PABM)")
  }
  
  n = nrow(A)
  
  if(model=="ER"){
    params = list(n=n, p = sum(A) / (n*(n-1)))
  }else if(model=="CL"){
    params = list(psi = colSums(A)/sqrt(sum(A)))
  }else if(model=="SBM"){
    K = length(unique(C))
    B <- matrix(0, nrow=K, ncol=K)
    
    G <- lapply(1:K, function(k) which(C == k))
    nk <- sapply(G, 'length')
    
    for(k in 1:K){
      for(l in k:K){
        B[k,l] <- B[l,k] <- sum(A[G[[k]], G[[l]]])/(nk[k]*nk[l])
      }
    }
    
    diag(B) <- diag(B)*nk/(nk-1)
    B[!is.finite(B)] <- 1e-6
    params = list(B=B, C=C)
  }else if(model=="DCBM"){
    out = fast.DCBM.est(A, C) # from NETCROP paper
    psi = out$psi
    B = out$Bsum
    params = list(psi=psi, B=B, C=C)
  }else if(model=="PABM"){
    foo = f.PA(A,C)
    M = foo$M
    O = foo$O
    lambda <- matrix(NA, nrow=n, ncol = max(C))
    for (i in 1:n){
      s <- C[i]
      for (r in 1:max(C)){
        lambda[i,r] <- M[i,r]/sqrt(O[s,r])} 
    }
    params = list(C=C, lambda = lambda)
  }
  
  return(params)
}

#' @title Split network with edge sampling
#' @description Splits the inputted adjacency matrix into two independent parts
#' @param A adjacency matrix
#' @param gamma parameter in data-splitting
#' @return Independent "split" adjacency matrices
#' @export
edge_sample <- function(A, gamma=0.5){
  n = dim(A)[1]
  Z <- matrix(rbinom(n^2, A, gamma), nrow=n)
  
  Z[lower.tri(Z)] <- 0
  Z = Z + t(Z)
  diag(Z) <- 0
  
  return(list(Z = Z, Y = A-Z))
}

#' @title Universal Inference e-value for a single run
#' @description Computes Universal Inference e-value on split network
#' @param A adjacency matrix
#' @param h0 null model
#' @param h1 alternative model
#' @param h0K number of communities in null hypothesis (ignored for ER or CL models)
#' @param h1K number of communities in alternative hypothesis (ignored for CL model)
#' @param gamma parameter in edge-splitting
#' @return e-value
#' @export
eval <- function(A, h0=c("ER", "CL", "SBM", "DCBM", "PABM"), h1=c("CL", "SBM", "DCBM", "PABM"), h0K, h1K, gamma=0.5){
  
  if(h0==h1 && h0K==h1K){
    stop("Null and alterantive hypotheses must be differnet.")
  }
  
  if(h1=="ER"){
    stop("Alternative hypothesis cannot be Erdos-Renyi.")
  }
  
  
  if(!(h0 %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
    Popularity adjusted block model (=PABM)")
  }
  
  if(!(h1 %in% c("CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
         Popularity adjusted block model (=PABM)")
  }
  
  out = edge_sample(A, gamma)
  Y = out$Y
  Z = out$Z
  
  if(h1=="CL"){
    paramsY = estparam(Y, "CL")
    paramsY$psi = paramsY$psi * sqrt(gamma)/ sqrt(1-gamma)
  }else if(h1=="SBM" || h1=="DCBM"){
    
    if(h1 == "SBM")  CY = randnet::reg.SP(A = Y, K = h1K, tau = 1)$cluster  # regular spectral clustering
    if(h1 == "DCBM") CY = randnet::reg.SSP(A = Y, K = h1K, tau = 1)$cluster # spherical spectral clustering
    
    paramsY = estparam(Y, h1, CY)
    paramsY$B = paramsY$B * gamma / (1-gamma)
  }else if(h1=="PABM"){
    CY = randnet::reg.SSP(A = Y, K = h1K, tau = 1)$cluster
    paramsY = estparam(Y, "PABM", CY)
    paramsY$lambda = paramsY$lambda * sqrt(gamma) / sqrt((1-gamma))
  }
  
  # Evaluate log-likelihood using parameters from Y on network Z
  L1 = llike(Z, h1, paramsY)
  
  # Estimate parameters under H0
  if(h0=="ER"){
    paramsZ = estparam(Z, "ER")
  }else if(h0=="CL"){
    paramsZ = estparam(Z, "CL")
  }else if(h0=="SBM"){
    CZ = randnet::reg.SP(A = Z, K = h0K, tau = 1)$cluster # regular spectral clustering
    paramsZ = estparam(Z, h0, CZ)
  }else if(h0=="DCBM"){
    CZ = randnet::reg.SSP(A = Z, K = h0K, tau = 1)$cluster # spherical spectral clustering
    paramsZ = estparam(Z, h0, CZ)
  }else if(h0=="PABM"){
    CZ = randnet::reg.SSP(A = Z, K = h0K, tau = 1)$cluster # spherical spectral clustering
    paramsZ = estparam(Z, h0, CZ)
  }
  # Evaluate log-likelihood using Z and parameters from Z
  L0 = llike(Z, h0, paramsZ)
  
  # Return the quotient, but upper bound it at ub for numerical stability
  ub = 1e100
  if((L1-L0) > log(ub)){
    return(ub)
  }else{
    return(exp(L1 - L0))
  }
  
}


#' @title Universal Inference e-value for a multiple runs
#' @description Computes Universal Inference e-value on split network
#' @param A adjacency matrix
#' @param h0 null model
#' @param h1 alternative model
#' @param h0K number of communities in null hypothesis (ignored for ER or CL models)
#' @param h1K number of communities in alternative hypothesis (ignored for CL model)
#' @param gamma parameter in edge-splitting
#' @param nreps number of repetitions of the data split
#' @param ncores number of cores for parallel computing
#' @return Mean e-value of nreps data splits
#' @export
eval_mc <- function(A, h0=c("ER", "CL", "SBM", "DCBM", "PABM"), h1=c("CL", "SBM", "DCBM", "PABM"), 
                 h0K, h1K, gamma=0.5, nreps, ncores){
  
  apply_fun <- function(i){
    
    err <- try(eval(A, h0, h1, h0K, h1K, gamma),TRUE)
    
    if(class(err)=="try-error"){
      return(NA)
    }else{
      return(err)
    }
    
  }
  
  out <- unlist(parallel::mclapply(1:nreps, apply_fun, mc.cores = ncores))
  
  if(mean(is.na(out)) > 0.10){
    stop("More than 10% of the iterations yielded an error.")
  }else{
    # Return the mean (and remove missing values)
    return(mean(out, na.rm = TRUE))
  }
  
}




eval_log <- function(A, h0=c("ER", "CL", "SBM", "DCBM", "PABM"), h1=c("CL", "SBM", "DCBM", "PABM"), h0K, h1K, gamma=0.5){
  
  if(h0==h1 && h0K==h1K){
    stop("Null and alterantive hypotheses must be differnet.")
  }
  
  if(h1=="ER"){
    stop("Alternative hypothesis cannot be Erdos-Renyi.")
  }
  
  
  if(!(h0 %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
    Popularity adjusted block model (=PABM)")
  }
  
  if(!(h1 %in% c("CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
         Popularity adjusted block model (=PABM)")
  }
  
  out = edge_sample(A, gamma)
  Y = out$Y
  Z = out$Z
  
  if(h1=="CL"){
    paramsY = estparam(Y, "CL")
    paramsY$psi = paramsY$psi * sqrt(gamma)/ sqrt(1-gamma)
  }else if(h1=="SBM" || h1=="DCBM"){
    
    if(h1 == "SBM")  CY = randnet::reg.SP(A = Y, K = h1K, tau = 1)$cluster  # regular spectral clustering
    if(h1 == "DCBM") CY = randnet::reg.SSP(A = Y, K = h1K, tau = 1)$cluster # spherical spectral clustering
    
    paramsY = estparam(Y, h1, CY)
    paramsY$B = paramsY$B * gamma / (1-gamma)
  }else if(h1=="PABM"){
    CY = randnet::reg.SSP(A = Y, K = h1K, tau = 1)$cluster
    paramsY = estparam(Y, "PABM", CY)
    paramsY$lambda = paramsY$lambda * sqrt(gamma) / sqrt((1-gamma))
  }
  
  # Evaluate log-likelihood using parameters from Y on network Z
  L1 = llike(Z, h1, paramsY)
  
  # Estimate parameters under H0
  if(h0=="ER"){
    paramsZ = estparam(Z, "ER")
  }else if(h0=="CL"){
    paramsZ = estparam(Z, "CL")
  }else if(h0=="SBM"){
    CZ = randnet::reg.SP(A = Z, K = h0K, tau = 1)$cluster # regular spectral clustering
    paramsZ = estparam(Z, h0, CZ)
  }else if(h0=="DCBM"){
    CZ = randnet::reg.SSP(A = Z, K = h0K, tau = 1)$cluster # spherical spectral clustering
    paramsZ = estparam(Z, h0, CZ)
  }else if(h0=="PABM"){
    CZ = randnet::reg.SSP(A = Z, K = h0K, tau = 1)$cluster # spherical spectral clustering
    paramsZ = estparam(Z, h0, CZ)
  }
  # Evaluate log-likelihood using Z and parameters from Z
  L0 = llike(Z, h0, paramsZ)
  
  # Return the quotient, but upper bound it at ub for numerical stability
  return(L1 - L0)
}

eval_mc_log <- function(A, h0=c("ER", "CL", "SBM", "DCBM", "PABM"), h1=c("CL", "SBM", "DCBM", "PABM"), 
                    h0K, h1K, gamma=0.5, nreps, ncores){
  
  apply_fun <- function(i){
    
    err <- try(eval_log(A, h0, h1, h0K, h1K, gamma),TRUE)
    
    if(class(err)=="try-error"){
      return(NA)
    }else{
      return(err)
    }
    
  }
  
  out <- unlist(parallel::mclapply(1:nreps, apply_fun, mc.cores = ncores))
  
  if(mean(is.na(out)) > 0.10){
    stop("More than 10% of the iterations yielded an error.")
  }else{
    # Return the median (and remove missing values)
    return(median(out, na.rm = TRUE))
  }
  
}


