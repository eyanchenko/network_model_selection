## Paper: Universal inference for model selection on networks
## Author: Eric Yanchenko
## University: Akita International University
## Date: June 29, 2026


source("netcrop.R")
source("functions.R")

library(ggplot2)
library(dplyr)
library(ggpubr)

## Setting 2 (a) - H0: DCBM w/ K=1 (CL) vs. H1: DCBM w/ K=2
## Increasing community structure strength

niter = 200
nreps = 100


n = 1000
alpha = 0.25
beta.seq = seq(0, 0.95, length=21)
delta=0.25
K = 2

methods = c("e-value0.4", "e-value0.5", "e-value0.6", "NETCROP", "ECV")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(beta.seq)), 
             Method=rep(methods, niter*length(beta.seq)), 
             beta=0, rej=0, time=0)

cnt = 1
for(beta in beta.seq){
  
  B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))
  pi = c(delta, 1-delta)
  
  for(iter in 1:niter){
    
    df$beta[cnt:(cnt+length(methods)-1)] <- beta
    
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.25, 0.75)
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    # E-value with gamma = 0.4
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "CL", "DCBM", 2, 2, 0.4, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "CL", "DCBM", 2, 2, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej
    
    # E-value with gamma = 0.6
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "CL", "DCBM", 2, 2, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses DCBM with 2 community
    df$time[cnt+3] <- system.time(rej <- as.numeric(netcrop_dcbm_sims(A, c(1,2), s=5, o=100, R=1, loss="l2")=="DCSBM-2"))[3]
    df$rej[cnt+3] <- rej
    
    # Reject if it chooses DCBM with 2 community (ECV)
    df$time[cnt+4] <- system.time(rej <- as.numeric(ecv_dcbm_sims(A, c(1,2))=="DCSBM-2"))[3]
    df$rej[cnt+4] <- rej
    
    cnt = cnt + length(methods)
    
    save(df, file="df_dcbm_K2_beta_062626.RData")
    
  }
  print(beta)
}

load("df_dcbm_K2_beta_062626.RData")

df_plot <- df %>% group_by(Method, beta) %>% summarize(rej = mean(rej, na.rm=TRUE), time=mean(time))

p1 <- ggplot(df_plot, aes(x=beta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(beta))+
  ylab("Rejection Rate")+
  theme_bw()+
  theme(text = element_text(size = 16))
p1


## Setting 2 (b) - H0: DCBM w/ K=1 (CL) vs. H1: DCBM w/ K=2
## Increasing community size

n = 1000
alpha = 0.25
beta = 0.7
delta.seq = seq(0.01, 0.5, length=20)
K = 2

B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))

methods = c("e-value0.4", "e-value0.5", "e-value0.6", "NETCROP", "ECV")


df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(delta.seq)), 
             Method=rep(methods, niter*length(delta.seq)), 
             delta=0, rej=0, time=0)

cnt = 1


for(delta in delta.seq){
  
  pi = c(delta, 1-delta)
  
  for(iter in 1:niter){
    
    df$delta[cnt:(cnt+length(methods)-1)] <- delta
  
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.25, 0.75)
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    # E-value with gamma = 0.1
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "CL", "DCBM", 2, 2, 0.4, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "CL", "DCBM", 2, 2, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej
    
    # E-value with gamma = 0.9
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "CL", "DCBM", 2, 2, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses SBM with 2 community
    df$time[cnt+3] <- system.time(
      rej <- as.numeric(netcrop_dcbm_sims(A, c(1,2), s=5, o=100, R=1, loss="l2")=="DCSBM-2")
    )[3]
    df$rej[cnt+3] <- rej
    
    # Reject if it chooses DCBM with 2 community (ECV)
    df$time[cnt+4] <- system.time(
      rej <- as.numeric(ecv_dcbm_sims(A, c(1,2))=="DCSBM-2")
    )[3]
    df$rej[cnt+4] <- rej
    
    cnt = cnt + length(methods)
    save(df, file="df_dcbm_K2_delta_062626.RData")
    
  }
  print(delta)
}

load("df_dcbm_K2_delta_062626.RData")

df_plot <- df %>% group_by(Method, delta) %>% summarize(rej = mean(rej, na.rm=TRUE), time=mean(time))

p2 <- ggplot(df_plot, aes(x=delta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("")+
  theme_bw()+
  theme(text = element_text(size = 16))
p2

ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")
ggsave("dcbm_K2_062626.pdf", height=4, width=8, unit="in")


## Setting 2 (c) - H0: DCBM w/ K=4 vs. H1: DCBM w/ K=5
## Increasing community size

n = 1000
alpha = 0.90
beta = 0.90
delta.seq = seq(0.01, 0.20, length=20)
K = 5

B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))

methods = c("e-value0.4", "e-value0.5", "e-value0.6", "NETCROP", "ECV")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(delta.seq)), 
             Method=rep(methods, niter*length(delta.seq)), 
             delta=0, rej=0, time=0)

cnt = 1


for(delta in delta.seq){
  
  pi = c(rep(1/(K-1) - delta/(K-1), K-1), delta)

  for(iter in 1:niter){
    
    df$delta[cnt:(cnt+length(methods)-1)] <- delta
    
    C  = sample(1:K, n, TRUE, pi)
    psi <- runif(n, 0.25, 0.75)
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    # E-value with gamma = 0.4
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "DCBM", K-1, K, 0.4, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "DCBM", K-1, K, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej
    
    # E-value with gamma = 0.6
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "DCBM", "DCBM", K-1, K, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses DCBM with 5 communities
    df$time[cnt+3] <- system.time(rej <- as.numeric(netcrop_dcbm_sims(A, c(4,5), s=5, o=100, R=1, loss="l2")=="DCSBM-5"))[3]
    df$rej[cnt+3] <- rej
    
    # Reject if it chooses DCBM with 5 communities (ECV)
    df$time[cnt+4] <- system.time(
      rej <- as.numeric(ecv_dcbm_sims(A, c(4,5))=="DCSBM-5")
    )[3]
    df$rej[cnt+4] <- rej
    

    cnt = cnt + length(methods)
    save(df, file="df_dcbm_K5_delta_062626.RData")
    
  }
  print(delta)
}

load("df_dcbm_K5_delta_062626.RData")

df_plot <- df %>% group_by(Method, delta) %>% summarize(rej = mean(rej), time=mean(time))

p1 <- ggplot(df_plot, aes(x=delta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Rejection Rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 16))
p1
ggsave("dcbm_K5_delta_062626.pdf", height=4, width=6, unit="in")


