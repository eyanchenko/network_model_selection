## Data fission and UI for model selection on networks
## Eric Yanchenko
## May 1, 2025


source("~/Documents/Research/network_model_selection/functions.R")

library(ggplot2)
library(dplyr)


## Setting 4 (a) -  H0: DCBM w/ K=2  vs. H1: PABM w/ K=2
niter = 200
nreps  = 100
n = 1000
beta = 3
delta = 0.25
K = 2
C   = c(rep(1, delta*n), rep(2, (1-delta)*n))
cat = c(rep(1,delta*n/2), rep(2,delta*n/2), rep(1,(1-delta)*n/2), rep(2,(1-delta)*n/2))

nu.seq = seq(0.10, 0.25, length=16)

methods = c("e-value0.2", "e-value0.4", "e-value0.5", "e-value0.6", "e-value0.8")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(nu.seq)), 
             Method=rep(methods, niter*length(nu.seq)), 
             nu=0, rej=0, time=0)

cnt = 1
for(nu in nu.seq){
  a = c(0.5+nu/2, 0.5-nu/2)
  b = c(0.5-nu/2, 0.5+nu/2)
  
  lambda <- matrix(0, n, K)
  
  for(i in 1:n){
    for(r in 1:K){
      
      if(C[i]==r){
        val = a[cat[i]] / sqrt(beta/(1+beta))
      }else{
        val = b[cat[i]] / sqrt(1/(1+beta))
      }
      
      lambda[i,r] <- val
    }
  }
  
  params = list(lambda=lambda, C=C)
  
  for(iter in 1:niter){
    
    df$nu[cnt:(cnt+length(methods)-1)] <- nu
    A <- generateA("PABM", params)
    
    # E-value with gamma = 0.2
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "PABM", K, K, 0.2, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.4
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "PABM", K, K, 0.4, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej

    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "PABM", K, K, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # E-value with gamma = 0.6
    # Reject if e-value is greater than 20.
    df$time[cnt+3] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "PABM", K, K, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+3]  <- rej
    
    # E-value with gamma = 0.8
    # Reject if e-value is greater than 20.
    df$time[cnt+4] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "PABM", K, K, 0.8, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+4]  <- rej
    
    cnt = cnt + length(methods)
    
    save(df, file="~/Documents/Research/network_model_selection/Results/df_dcbm_pabm_K2_nu.RData")
  }
  print(nu)
}


load("~/Documents/Research/network_model_selection/Results/df_dcbm_pabm_K2_nu.RData")

df_plot <- df %>% group_by(Method, nu) %>% summarize(rej = mean(rej, na.rm=TRUE), time=mean(time))

p1 <- ggplot(df_plot, aes(x=nu, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Rejection Rate")+
  theme_bw()
p1

ggsave("~/Documents/Research/network_model_selection/Figures/dcbm_pabm_K5_nu.pdf", height=4, width=6, unit="in")



