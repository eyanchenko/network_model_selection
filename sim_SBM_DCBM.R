## Data fission and UI for model selection on networks
## Eric Yanchenko
## May 1, 2025


source("~/Documents/Research/network_model_selection/netcrop.R")
source("~/Documents/Research/network_model_selection/functions.R")

library(ggplot2)
library(dplyr)


## Setting 3 (a) -  H0: SBM w/ K=2  vs. H1: DCBM w/ K=2
## Increasing variance of degree parameters
niter = 100
n = 1000
alpha = 0.5
beta = 0.20
pp = 0.5
prop = 0.75
K = 2
B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))
pi = c(prop, 1-prop)

nu.seq <- seq(0.25, 0.5, length=11)

df <- tibble(iter = rep(rep(1:niter,each=2), length(nu.seq)), 
             Method=rep(c("e-value", "NETCROP"), niter*length(nu.seq)), 
             nu=0, rej=0, time=0)

cnt = 1
for(nu in nu.seq){
  
  for(iter in 1:niter){
    
    df$nu[c(cnt,cnt+1)] <- nu
    
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.5 - nu/2, 0.5 + nu/2)
    
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    
    # Reject if e-value is greater than 10.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, pp, nreps = 100, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # Reject if it chooses DCBM with 2 community
    df$time[cnt+1] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=K, s=5, o=100, R=1, loss="l2")$l2.model=="DCSBM-2")
    )[3]
    df$rej[cnt+1] <- rej
    
    cnt = cnt + 2
    
    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K2_nu.RData")
    
  }
  print(nu)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K2_nu.RData")

df_plot <- df %>% group_by(Method, nu) %>% summarize(rej = mean(rej), time=mean(time))

p1 <- ggplot(df_plot, aes(x=nu, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Rejection Rate")+
  theme_bw()
p1

ggsave("~/Documents/Research/network_model_selection/Figures/sbm_dcbm_K2_nu.pdf", height=4, width=6, unit="in")


p2 <- ggplot(df_plot, aes(x=nu, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(nu))+
  ylab("Time (sec)")+
  theme_bw()
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")


## Setting 3 (b) -  H0: SBM w/ K=5  vs. H1: DCBM w/ K=5
## Increasing variance of degree parameters
niter = 100
n = 1000
alpha = 0.5
beta = 0.20
pp = 0.5
K = 5
B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))
pi = c(0.40, 0.15, 0.15, 0.15, 0.15)

nu.seq <- seq(0.25, 0.5, length=11)

df <- tibble(iter = rep(rep(1:niter,each=2), length(nu.seq)), 
             Method=rep(c("e-value", "NETCROP"), niter*length(nu.seq)), 
             nu=0, rej=0, time=0)

cnt = 1
for(nu in nu.seq){
  
  for(iter in 1:niter){
    
    df$nu[c(cnt,cnt+1)] <- nu
    
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.5 - nu/2, 0.5 + nu/2)
    
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    
    # Reject if e-value is greater than 10.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "DCBM", K, K, pp, nreps = 100, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # Reject if it chooses DCBM with 2 community
    df$time[cnt+1] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=K, s=5, o=100, R=1, loss="l2")$l2.model=="DCSBM-5")
    )[3]
    df$rej[cnt+1] <- rej
    
    cnt = cnt + 2
    
    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K5_nu.RData")
    
  }
  print(nu)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_dcbm_K5_nu.RData")

df_plot <- df %>% group_by(Method, nu) %>% summarize(rej = mean(rej), time=mean(time))

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
