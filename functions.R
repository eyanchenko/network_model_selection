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
generateA <- function(model=c("ER", "CL", "SBM", "DCBM", "RDPG", "MMSBM", "LSM"), params){
  
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "RDPG", "MMSBM", "LSM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
         Random Dot Product Graph (RDPG), mixed-membership SBM (MMSBM),
         Latent space model (LSM)")
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
  }else if(model=="RDPG"){
    n = params$n
    d = params$d
    xi = params$xi
    
    X = matrix(runif(n*d), nrow=n, ncol=d)
    XXt = X%*%t(X)
    P = xi*XXt / max(XXt)
  }else if(model=="MMSBM"){
    PI = params$PI
    B  = params$B
    n  = nrow(PI)
    K  = ncol(PI)
    
    P = matrix(0, nrow = n, ncol = n)
    
    for(i in 2:n){
      for(j in 1:(i-1)){
        ci = sample(1:K, 1, FALSE, PI[i, ])
        cj = sample(1:K, 1, FALSE, PI[j, ])
        
        P[i,j] = P[j,i] = B[ci, cj]
      }
    }
  } else if (model=="LSM"){
    n     = params$n
    d     = params$d
    alpha = params$alpha
  
    Z = matrix(rnorm(n*d), nrow=n, ncol=d)
    
    P = matrix(0, nrow = n, ncol = n)
    
    for(i in 2:n){
      for(j in 1:(i-1)){
        eta = alpha - sum((Z[i, ] - Z[j, ])^2)
        P[i,j] = P[j,i] = 1/(1+exp(-eta))
      }
    }
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
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM", "RDPG", "MMSBM", "LSM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), Popularity adjusted block model (PABM), degree-corrected block model (DCBM),
         Random Dot Product Graph (RDPG)", "Mixed-membership SBM (MMSBM), 
         Latent-space model (LSM)")
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
    C = params$C
    lambda = params$lambda
    n = length(C)
    
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
     
  }else if (model=="RDPG"){
    X = params$X
    P = X%*%t(X)  
    
    eps = 1e-6 # tolerance
    P[P < eps] <- eps
    P[P > 1-eps] <- 1 - eps
    
  }else if (model=="MMSBM"){
    PI = params$PI
    B  = params$B
    n  = nrow(PI)
    K  = ncol(PI)
    P  = matrix(0, nrow = n, ncol = n)
    for(i in 2:n){
      for(j in 1:(i-1)){
        P[i,j] = P[j,i] = t(PI[i,])%*%B%*%PI[j,]
      }
    }
  }else if (model=="LSM"){
    alpha = params$alpha
    Z = params$Z
    
    P  = matrix(0, nrow = n, ncol = n)
    for(i in 2:n){
      for(j in 1:(i-1)){
        eta = alpha - sum((Z[i, ] - Z[j, ])^2)
        P[i,j] = P[j,i] = 1/(1+exp(-eta))
      }
    }
  }
  
  eps = 1e-6
  P[P < eps] <- eps
  P[P > 1-eps] <- 1-eps
  
  return(
    -sum(P[upper.tri(P)]) + sum(A[upper.tri(A)] * log(P[upper.tri(P)]))
  )
}

# Needed for PABM estimation
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
#' @param d dimension of RDPG latent space(if applicable)
#' @param mmK number of communities for MMSBM (if applicable)
#' @return estimated parameter values
#' @export
estparam <- function(A, model, C, d, mmK){
  if(!(model %in% c("ER", "CL", "SBM", "DCBM", "PABM", "RDPG", "MMSBM", "LSM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
    Popularity adjusted block model (=PABM), Random Dot Product Graph (RDPG), mixed-membership SBM (MMSBM), Latent-space model (LSM)")
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
  }else if (model=="RDPG"){
    #deg = rowSums(A) # if you want to use the Laplacian
    #D <- sparseMatrix(i = 1:n, j = 1:n, x = 1/sqrt(deg))
    #L <- tcrossprod(crossprod(D, A), D)
    
    eig <- irlba::irlba(A, nv = d)
    U <- eig$v
    sigma.half <- diag(sqrt(abs(eig$d)), nrow = d, ncol = d)
    X <- tcrossprod(U, sigma.half)
    params = list(X=X)
  }else if(model=="MMSBM"){
    # Convert to edge list
    edge_list = tibble(v1=numeric(n*(n-1)/2), v2=0, Y=0)
    idx = 1
    for(i in 2:n){
      for(j in 1:(i-1)){
        edge_list[idx, 1] = i
        edge_list[idx, 2] = j
        edge_list[idx, 3] = A[i,j]
        idx = idx + 1
      }
    }
    
    out <- mmsbm(Y~1, data.dyad = edge_list, senderID = "v2", receiverID = "v1",
                 n.blocks = mmK, directed = FALSE)
    
    PI = as.matrix(t(out$MixedMembership))
    hold = capture.output(outt <- summary(out))
    B  = outt$`Blockmodel Matrix`
    
    params = list(PI=PI, B=B)
  }else if (model=="LSM"){
    out <- LSM.PGD(A, d, step.size=0.3,niter=50,trace=0)
    params = list(alpha = mean(out$alpha), Z = out$Z)
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

#' @title Universal Inference e-value for a single run
#' @description Computes Universal Inference e-value on split network
#' @param A adjacency matrix
#' @param h0 null model
#' @param h1 alternative model
#' @param h0K number of communities (or dim of latent space) in null hypothesis
#' @param h1K number of communities (or dim of latent space) in alternative hypothesis
#' @param theta parameter in data-splitting
#' @return e-value
#' @export
eval <- function(A, h0=c("ER", "CL", "SBM", "DCBM", "PABM", "RDPG"), h1=c("CL", "SBM", "DCBM", "RDPG", "MMSBM", "LSM"), 
                 h0K, h1K, theta=0.5){
  
  if(h0==h1 && h0K==h1K){
    stop("Null and alterantive hypotheses must be differnet.")
  }
  
  if(h1=="ER"){
    stop("Alternative hypothesis cannot be Erdos-Renyi.")
  }
  
  
  if(!(h0 %in% c("ER", "CL", "SBM", "DCBM", "PABM", "RDPG", "LSM"))){
    stop("Please enter a valid model name. Erdos-Renyi (=ER), Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
    Popularity adjusted block model (=PABM),
         random dot product graph (RDPG), Latent-space model (LSM)")
  }
  
  if(!(h1 %in% c("CL", "SBM", "DCBM", "PABM", "RDPG", "MMSBM", "LSM"))){
    stop("Please enter a valid model name. Chung-Lu (=CL),
    Stochastic-block model (=SBM), degree-corrected block model (DCBM),
         Popularity adjusted block model (=PABM),
         random dot product graph (RDPG), mixed-membership SBM (MMSBM),
         Latent-space model (LSM)")
  }
  
  out = fission(A, theta)
  Y = out$Y
  Z = out$Z
  
  if(h1=="CL"){
    paramsY = estparam(Y, "CL")
    paramsY$psi = paramsY$psi * sqrt(theta)/ sqrt(1-theta)
  }else if(h1=="SBM" || h1=="DCBM"){
    
    if(h1 == "SBM")  CY = randnet::reg.SP(A = Y, K = h1K, tau = 1)$cluster  # regular spectral clustering
    if(h1 == "DCBM") CY = randnet::reg.SSP(A = Y, K = h1K, tau = 1)$cluster # spherical spectral clustering
    
    paramsY = estparam(Y, h1, CY)
    paramsY$B = paramsY$B * theta / (1-theta)
  }else if(h1=="PABM"){
    CY = randnet::reg.SSP(A = Y, K = h1K, tau = 1)$cluster
    paramsY = estparam(Y, "PABM", CY)
    paramsY$lambda = paramsY$lambda * sqrt(theta) / sqrt((1-theta))
  }else if(h1=="RDPG"){
    XY = estparam(Y, "RDPG", NULL, h1K)$X * sqrt(theta) / sqrt((1-theta))
    paramsY = list(X=XY)
  }else if(h1 == "MMSBM"){
    paramsY = estparam(A, "MMSBM", mmK = h1K)
    paramsY$B = paramsY$B * theta/(1-theta)
  }else if(h1 == "LSM"){
    paramsY = estparam(A, "LSM", d=h1K) # Does this need to be changed if theta != 0.5?
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
  }else if(h0=="RDPG"){
    XZ = estparam(Z, "RDPG", d=h0K)$X
    paramsZ = list(X=XZ)
  }else if(h0=="LSM"){
    paramsZ = estparam(Z, "LSM", d=h0K)
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
#' @param h0K number of communities (or dim of latent space) in null hypothesis
#' @param h1K number of communities (or dim of latent space) in alternative hypothesis
#' @param theta parameter in data-splitting
#' @param nreps number of repetitions of the data split
#' @param ncores number of cores for parallel computing
#' @return Median e-value of nreps data splits
#' @export
eval_mc <- function(A, h0=c("ER", "CL", "SBM", "DCBM", "RDPG"), h1=c("CL", "SBM", "DCBM", "RDPG", "MMSBM", "LSM"), 
                 h0K, h1K, theta=0.5, nreps, ncores){
  
  apply_fun <- function(i){
    
    err <- try(eval(A, h0, h1, h0K, h1K, theta),TRUE)
    
    if(class(err)=="try-error"){
      return(NA)
    }else{
      return(err)
    }
    
  }
  
  out <- unlist(parallel::mclapply(1:nreps, apply_fun, mc.cores = ncores))
  
  if(mean(is.na(out)) > 0.05){
    stop("More than 5% of the iterations yielded an error.")
  }else{
    # Return the median (and remove missing values)
    return(median(out, na.rm = TRUE))
  }
  
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
    A.i <- generateA("ER", list(n=n, p=p.hat))
    A.i.prime <- (A.i-P.hat)/sqrt((n-1)*p.hat*(1-p.hat))
    princ.eigen.i <- RSpectra::eigs_sym(A.i.prime,1,which="LA")[[1]]
    emp.stats[i] <- n^(2/3)*(princ.eigen.i-2)
  }
  
  mu.theta <- mean(emp.stats)
  sigma.theta <- sqrt(var(emp.stats))
  
  theta.prime <- mu.tw + (obs.stat-mu.theta)/sigma.theta * sigma.tw
  return(RMTstat::ptw(theta.prime, beta=1, lower.tail = FALSE))
  
}


## Yanchenko and Sengupta (2024) bootstrap method

# Greedy community detection algorithm
greedy <- function(A, K, runs=2, max.iter=100){
  
  n = dim(A)[1]
  
  # Find max e2d2 and corresponding labels
  apply.fun <- function(i){
    # Initialize communities
    comm = sample(1:K,n,replace=T)
    comm_old = comm
    
    # Community sizes
    nc = numeric(K)
    for(i in 1:K){nc[i] <- sum(comm==i)}
    
    
    # Initial values
    M   = sum(A)/2
    Xin = 0
    for(i in 1:K){Xin = Xin + sum(A[comm==i, comm==i])/2}
    Xout = M - Xin
    
    m   = 0.5*n*(n-1)
    m_in = 0
    for(i in 1:K){m_in = m_in + 0.5*nc[i]*(nc[i]-1)}
    m_out = m - m_in
    
    e2d2 = Xin/m_in - Xout/m_out
    
    run = T
    iter = 1
    while(run && iter < max.iter){
      run = F
      # Randomize node order
      node.order = sample(1:n)
      for(i in node.order){
        
        #Find community neighbors of node i
        neighs = sample(unique(comm[as.logical(A[i,])]))
        e2d2_hold = numeric(length(neighs))
        
        # Compute change in e2d2 for all neighboring communities
        Xlost = sum(A[i, comm==comm[i]])
        m_in_lost = nc[comm[i]] - 1
        
        apply.fun2 <- function(j){
          # Swap community of node i
          Xin_i = Xin - Xlost + sum(A[i, comm==neighs[j]])
          Xout_i = M - Xin_i
          
          m_in_i = m_in - m_in_lost + nc[neighs[j]] - 1*(neighs[j]==comm[i])
          m_out_i = m - m_in_i
          
          return(Xin_i/m_in_i - Xout_i/m_out_i)
        }
        
        e2d2_hold <- unlist(lapply(1:length(neighs), apply.fun2))
        
        # for(j in 1:length(neighs)){
        #   # Swap community of node i
        #   Xin_i = Xin - Xlost + sum(A[i, comm==neighs[j]])
        #   Xout_i = M - Xin_i
        # 
        #   m_in_i = m_in - m_in_lost + nc[neighs[j]] - 1*(neighs[j]==comm[i])
        #   m_out_i = m - m_in_i
        # 
        #   e2d2_hold[j] <- Xin_i/m_in_i - Xout_i/m_out_i
        # 
        # }
        
        # Swap node to whichever community yields largest e2d2 and update values
        ck = neighs[which.max(e2d2_hold)]
        
        if(length(ck)==0){
          ck = comm[i]
        }
        
        Xin = Xin - Xlost + sum(A[i, comm==ck])
        Xout = M - Xin
        
        m_in = m_in - m_in_lost + nc[ck] - 1*(ck==comm[i])
        m_out = m - m_in
        
        e2d2 <-  Xin/m_in - Xout/m_out
        nc[comm[i]] = nc[comm[i]] - 1
        nc[ck]      = nc[ck] + 1
        comm[i] <- ck
        
        
        
      }
      
      
      
      # Stop if no labels updated
      if(sum(comm==comm_old) < n){
        run = T
        comm_old = comm
      }
      
      iter = iter + 1
      
    }
    
    return(c(comm, e2d2/(M/(n*(n-1)/2))/K ))
  }
  
  out <- lapply(1:runs, apply.fun)
  
  # Keep run with max e2d2
  e2d2_vec <- numeric(runs)
  for(i in 1:runs){
    e2d2_vec[i] <- out[[i]][n+1]
  }
  
  keeper = which.max(e2d2_vec)
  
  return(list(e2d2=out[[keeper]][n+1], comm=out[[keeper]][1:n]))
  
}

# Compute test statistic
test.stat <- function(A, C){
  n = dim(A)[1]
  
  K = length(unique(C))
  
  if (K == 1 || K == n) {
    
    return(-1)
    
  }
  
  
  P = choose(n,2)
  P.in <- 0
  for(i in 1:K){P.in = P.in + choose(length(C[[i]]),2)} 
  P.out = P-P.in
  
  m.in <- 0
  for(i in 1:K){m.in = m.in + sum(A[C[[i]],C[[i]]])/2}
  m.out <- sum(A)/2-m.in
  
  p.in <- m.in/P.in 
  p.out <- m.out/P.out
  
  
  p.bar <- sum(A)/(n*(n-1))
  return((p.in - p.out)/p.bar/K)
  
  
}

# Bootstrap p-value
boot.pval <- function(A, C, nsim, null=c("ER", "CL"), rn=1){
  #A: observed adjacency matrix
  #nsim: # iterations
  #Null: What null hypothesis do you want to compare against?
  #ER: Erdos-Renyi, CL: Chung-Lu
  #given community detection algorithm
  
  # Compute test statistic
  ts = test.stat(A, C)

  n = dim(A)[1]
  
  if(null=="ER"){
    p.hat <- sum(A)/(n*(n-1))
    apply.fun <- function(i, nl){
      A.hat <- generateA("ER", params=list(n=n, p=p.hat))
    }
        
  }else if(null=="CL"){
    fit <- eigs(A, 1, which="LM")
    theta.hat <- abs(as.numeric(fit$vectors * sqrt(fit$values)))
    theta.hat[theta.hat > 1] <- 1
    apply.fun <- function(i, nl){
        theta.boot <- sample(theta.hat, n, replace=T)
        A.hat <- generateA("ER", list(psi=theta.boot))
  }
    
    # Only keep largest connected component
    A.graph <- as.undirected(graph.adjacency(A.hat, mode="undirected"))
    comps <- components(A.graph)
    big_id <- which.max(comps$csize)
    v_ids <- V(A.graph)[comps$membership == big_id]
    A.graph <- induced_subgraph(A.graph, v_ids)
    A.hat <- as.matrix(as_adjacency_matrix(A.graph))   
    K <- length(cluster_fast_greedy(A.graph))
    
    greedy(A.hat, K, runs=rn)$e2d2
  }
  
  emp.samp <- unlist(mclapply(1:nsim, apply.fun, nl=null, mc.cores = 6))
  
  
  pval <- mean(emp.samp > ts)
  
  return(pval)
}

# Asymptotic cut-off
asymp.cut.agn <- function(A, C, K=2, g0=0){
  n = dim(A)[1]
  
  ts <- test.stat(A, C)
  
  p.hat <- sum(A)/(n*(n-1))
  eps = 0.0001
  
  cutoff = (g0 + (log(K)/n)^(1/2)/p.hat/K) * (1+eps)
  
  return(as.numeric(ts > cutoff))
  
}









