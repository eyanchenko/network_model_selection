## Paper: Universal inference for model selection on networks
## Author: Eric Yanchenko
## Akita International University


source("~/Documents/Research/network_model_selection/netcrop.R")
source("~/Documents/Research/network_model_selection/functions.R")

library(ggplot2)
library(dplyr)


## Setting 3 (a) -  H0: SBM w/ K=2  vs. H1: DCBM w/ K=2
## Increasing variance of degree parameters
niter = 200
nreps = 100

### Are these the settings that we used?
n = 1000
alpha = 0.25
beta = 0.50
delta = 0.25
K = 2
B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))
pi = c(delta, 1-delta)


nu.seq <- seq(0, 0.8, length=21)

methods = c("e-value0.4", "e-value0.5", "e-value0.6", "NETCROP", "ECV")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(nu.seq)), 
             Method=rep(methods, niter*length(nu.seq)), 
             nu=0, rej=0, time=0)

cnt = 1
for(nu in nu.seq){
  
  for(iter in 1:niter){
    
    df$nu[cnt:(cnt+length(methods)-1)] <- nu
    
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.5 - nu/2, 0.5 + nu/2)
    
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    # E-value with gamma = 0.4
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, 0.1, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej

    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej

    # E-value with gamma = 0.6
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses DCBM with K communities
    df$time[cnt+3] <- system.time(
      rej <- as.numeric(netcrop_sbm_dcbm_sims(A, 2, s=5, o=100, R=1, loss="l2")=="DCSBM-2")
    )[3]
    df$rej[cnt+3] <- rej
    
    # Reject if it chooses DCBM with 2 communities (ECV)
    df$time[cnt+4] <- system.time(rej <- as.numeric(ecv_sbm_dcbm_sims(A, 2)=="DCSBM-2"))[3]
    df$rej[cnt+4] <- rej
    
    
    cnt = cnt + length(methods)
    
    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K2_nu.RData")
    
  }
  print(nu)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K2_nu.RData")

df_plot <- df %>% group_by(Method, nu) %>% summarize(rej = mean(rej, na.rm=TRUE), time=mean(time))

p1 <- ggplot(df_plot, aes(x=nu, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Rejection Rate")+
  theme_bw()
p1

ggsave("~/Documents/Research/network_model_selection/Figures/sbm_dcbm_K2_nu_111725.pdf", height=4, width=6, unit="in")


p2 <- ggplot(df_plot, aes(x=nu, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Time (sec)")+
  theme_bw()
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")


## Setting 3 (b) -  H0: SBM w/ K=5  vs. H1: DCBM w/ K=5
## Increasing variance of degree parameters

n = 1000
alpha = 0.25
beta = 0.90
K = 5
B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))
#pi = c(0.40, 0.15, 0.15, 0.15, 0.15)
pi = rep(1/K, K)

nu.seq <- seq(0, 0.75, length=21)

methods = c("e-value0.4", "e-value0.5", "e-value0.6", "NETCROP", "ECV")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(nu.seq)), 
             Method=rep(methods, niter*length(nu.seq)), 
             nu=0, rej=0, time=0)

cnt = 1
for(nu in nu.seq){
  
  for(iter in 1:niter){
    
    df$nu[cnt:(cnt+length(methods)-1)] <- nu
    
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.5 - nu/2, 0.5 + nu/2)
    
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    eval_mc(A, "SBM", "DCBM", K, K, 0.4, nreps = nreps, ncores = detectCores()-1)
    
    
    # E-value with gamma = 0.4
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, 0.4, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej

    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej

    # E-value with gamma = 0.6
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses DCBM with K communities
    df$time[cnt+3] <- system.time(
      rej <- as.numeric(netcrop_sbm_dcbm_sims(A, 5, s=5, o=100, R=1, loss="l2")=="DCSBM-5")
    )[3]
    df$rej[cnt+3] <- rej
    
    # Reject if it chooses DCBM with 5 communities (ECV)
    df$time[cnt+4] <- system.time(
      rej <- as.numeric(ecv_sbm_dcbm_sims(A, 5)=="DCSBM-5")
    )[3]
    df$rej[cnt+4] <- rej
    
    cnt = cnt + length(methods)
    
    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K5_nu.RData")
    
  }
  print(nu)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K5_nu.RData")

df_plot <- df %>% group_by(Method, nu) %>% summarize(rej = mean(rej, na.rm=TRUE), time=mean(time))

p1 <- ggplot(df_plot, aes(x=nu, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Rejection Rate")+
  theme_bw()
p1

ggsave("~/Documents/Research/network_model_selection/Figures/sbm_dcbm_K5_nu.pdf", height=4, width=6, unit="in")


p2 <- ggplot(df_plot, aes(x=nu, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Time (sec)")+
  theme_bw()
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")
