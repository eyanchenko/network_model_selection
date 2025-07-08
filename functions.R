## Data fission and UI for model selection on networks
## Eric Yanchenko
## May 1, 2025

#' @title Generate networks from models
#' @description Generates network from variance 
#' random-graph models
#' @param model random-graph model
#' @param params (list of) parameters corresponding to the particular model
#' @return adjacency matrix
#' @export
generateA <- function(model=c("ER", "CL", "SBM", "DCBM", "PABM"), params){
  
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM), 
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
llike <- function(A, model, params){
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM), 
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
    
  }
  
  return(
    -sum(P[upper.tri(P)]) + sum(A[upper.tri(A)] * log(P[upper.tri(P)]))
  )
}

#' @title Spectral clustering
#' @description Estimate community labels using Spectral clustering
#' @param A adjacency matrix
#' @param K number of communities
#' @return community labels
#' @export
spectral <- function(A, K){
  n = nrow(A)
  dd = colSums(A)
  D <- Matrix::sparseMatrix( i = 1:n, j = 1:n,
                     x = 1/sqrt(dd))
  
  L <- tcrossprod(crossprod(D, A), D) # normalized graph Laplacian
  
  U = irlba::partial_eigen(L, n = K, 
                           symmetric = T)$vectors
  nc = U[,1]^2
  for(i in 2:K){nc = nc + U[,i]^2}
  nc = sqrt(nc)
  eV = U / nc # Normalize to have rows with norm 1
  
  C <- kmeans(eV, K)$cluster # k-means cluster
  return(C)
}

#' @title Estimate model parameters
#' @description Estimates model parameters for random graph model
#' @param A adjacency matrix
#' @param model random-graph model
#' @param C community labels (if applicable)
#' @return estimated parameter values
#' @export
estparam <- function(A, model, C){
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM), 
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
    out = fast.DCBM.est(A, C) # From NETCROP paper
    psi = out$psi
    B = out$Bsum
    params = list(psi=psi, B=B, C=C)
  }else if(model=="PABM"){
    
  }
  
  
  return(params)
}

#' @title Data fission of network
#' @description Splits the inputted adjacency matrix into two independent parts
#' @param A adjacency matrix
#' @param theta parameter in data-splitting
#' @return Independent "split" adjacency matrices
#' @export
fission <- function(A, theta=0.5){
  n = dim(A)[1]
  Z <- matrix(rbinom(n^2, A, theta), nrow=n)
  
  Z[lower.tri(Z)] <- 0
  Z = Z + t(Z)
  diag(Z) <- 0
  
  return(list(Z = Z, Y = A-Z))
}

#' @title Universal Inference e-value
#' @description Computes Universal Inference e-value on split network
#' @param A adjacency matrix
#' @param h0 null model
#' @param h1 alternative model
#' @param h0K number of communities in null hypothesis
#' @param h1K number of communities in alternative hypothesis
#' @param theta parameter in data-splitting
#' @return e-value
#' @export
eval <- function(A, h0=c("ER", "CL", "SBM", "DCBM"), h1=c("CL", "SBM", "DCBM", "PABM"), 
                 h0K, h1K, theta=0.5){
  
  if(h0==h1 && h0K==h1K){
    stop("Null and alterantive hypotheses must be differnet.")
  }
  
  if(h1=="ER"){
    stop("Alternative hypothesis cannot be Erdos-Renyi.")
  }
  
  if(h0=="PABM"){
    stop("Null hypothesis cannot be PABM.")
  }
  
  if(!(h0 %in% c("ER", "CL", "SBM", "DCBM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM)")
  }
  
  if(!(h1 %in% c("CL", "SBM", "DCBM", "PABM"))){
    stop("Please enter a valid model name. Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM), 
    popularity adjusted block model (=PABM)")
  }
  
  out = fission(A, theta)
  Y = out$Y
  Z = out$Z
  
  if(h1=="CL"){
    paramsY = estparam(Y, "CL")
    paramsY$psi = paramsY$psi * sqrt(theta)/ sqrt(1-theta)
  }else{
    CY = spectral(Y, h1K)
    paramsY = estparam(Y, h1, CY)
    paramsY$B = paramsY$B * theta / (1-theta)
  }

  
  # Evaluate log-likelihood using parameters from Y on network Z
  L1 = llike(Z, h1, paramsY)
  
  # Estimate parameters under H0
  if(h0=="ER"){
    paramsZ = estparam(Z, "ER")
  }else if(h1=="CL"){
    paramsZ = estparam(Z, "CL")
  }else{
    CZ = spectral(Z, h0K)
    paramsZ = estparam(Z, h0, CZ)
  }
  
  # Evaluate log-likelihood using Z and parameters from Z
  L0 = llike(Z, h0, paramsZ)
  
  
  # Return the quotient
  return(exp(L1 - L0))
}

#' @title Community detection p-value
#' @description Computes p-value for testing against ER null using Bickel and Sarkar (2016) method
#' @param A adjacency matrix
#' @return p-value
#' @export
spectral.pval <- function(A){
  
  n=dim(A)[1]
  
  p.hat <- sum(A)/(n*(n-1))
    
  P.hat <- p.hat - p.hat*diag(1,n)
  A.prime <- (A-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    
  princ.eigen <- RSpectra::eigs_sym(A.prime,1,which="LA")[[1]]
    
  obs.stat <- n^(2/3)*(princ.eigen-2)
  return(RMTstat::ptw(obs.stat, beta=1, lower.tail = FALSE))
}

#' @title Community detection p-value (with bootstrap)
#' @description Computes p-value for testing against ER null using Bickel and Sarkar (2016) bootstrap method
#' @param A adjacency matrix
#' @return p-value
#' @export
spectral.adj.pval <- function(A){
  
  n=dim(A)[1]

  p.hat <- sum(A)/(n*(n-1))
  
  P.hat <- p.hat - p.hat*diag(1,n)
  A.prime <- (A-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
  
  princ.eigen <- RSpectra::eigs_sym(A.prime,1,which="LA")[[1]]
  
  obs.stat <- n^(2/3)*(princ.eigen-2)
  
  mu.tw <- -1.2065335745820 #from wikipedia
  sigma.tw <- sqrt(1.607781034581) #from wikipedia
  
  emp.stats <- numeric(50)
  
  for(i in 1:50){
    A.i <- generateER(n, p.hat)
    A.i.prime <- (A.i-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    princ.eigen.i <- RSpectra::eigs_sym(A.i.prime,1,which="LA")[[1]]
    emp.stats[i] <- n^(2/3)*(princ.eigen.i-2)
  }
  
  mu.theta <- mean(emp.stats)
  sigma.theta <- sqrt(var(emp.stats))
  
  theta.prime <- mu.tw + (obs.stat-mu.theta)/sigma.theta * sigma.tw
  return(RMTstat::ptw(theta.prime, beta=1, lower.tail = FALSE))
  
}

















