## Data fission and UI for model selection on networks
## Eric Yanchenko
## May 1, 2025


source("~/Documents/Research/network_model_selection/netcrop.R")
source("~/Documents/Research/network_model_selection/functions.R")

library(ggplot2)
library(dplyr)

## Setting 1 (a) - H0: SBM w/ K=1 (ER) vs. H1: SBM w/ K=2
## Increasing community structure strength

# I think the clustering algorithms aren't working great...
# Could include an oracle method as comparison where we use the known C
# It's not our fault if the labels can't be estimated well...

niter = 100
n = 1000
alpha = 0.1
beta.seq = seq(0.2, 0.6, length=21)
delta=0.25
nreps = 100
K = 2

methods = c("e-value0.1", "e-value0.5", "e-value0.9", "NETCROP")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(beta.seq)), 
             Method=rep(methods, niter*length(beta.seq)), 
             beta=0, rej=0, time=0)

cnt = 1
for(beta in beta.seq){
  
  B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))
  
  pi = c(delta, 1-delta)
  C  = sample(1:K, n, TRUE, pi)
  
  for(iter in 1:niter){
    
    df$beta[cnt:(cnt+length(methods)-1)] <- beta
    
    params = list(B=B, C=C)
    A <- generateA("SBM", params)

    # E-value with gamma = 0.1
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "ER", "SBM", 2, 2, 0.1, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "ER", "SBM", 2, 2, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej
    
    # E-value with gamma = 0.9
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "ER", "SBM", 2, 2, 0.9, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej

    # Reject if it chooses SBM with 2 community
    df$time[cnt+3] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=2, s=5, o=100, R=1, loss="l2")$l2.model=="SBM-2")
    )[3]
    df$rej[cnt+3] <- rej

    cnt = cnt + length(methods)

    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_K2_beta.RData")

  }
  print(beta)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_K2_beta.RData")

df_plot <- df %>% group_by(Method, beta) %>% summarize(rej = mean(rej, na.rm=TRUE), time=mean(time))

p1 <- ggplot(df_plot, aes(x=beta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(beta))+
  ylab("Rejection Rate")+
  theme_bw()
p1

ggsave("~/Documents/Research/network_model_selection/Figures/sbm_K2_beta.pdf", height=4, width=6, unit="in")

p2 <- ggplot(df_plot, aes(x=beta, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(beta))+
  ylab("Time (sec)")+
  theme_bw()

library(ggpubr)
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")

## Setting 1 (b) - H0: SBM w/ K=1 (ER) vs. H1: SBM w/ K=2
## Increasing community size

niter = 100
n = 1000
alpha = 0.1
beta = 0.4
delta.seq = seq(0.05, 0.5, length=10)
K = 2

B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))

methods = c("e-value0.1", "e-value0.5", "e-value0.9", "NETCROP")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(delta.seq)), 
             Method=rep(methods, niter*length(delta.seq)), 
             delta=0, rej=0, time=0)

cnt = 1


for(delta in delta.seq){
  
  pi = c(delta, 1-delta)
  C  = sample(1:K, n, TRUE, pi)
  
  for(iter in 1:niter){
    
    df$delta[cnt:(cnt+length(methods)-1)] <- delta
    
    params = list(B=B, C=C)
    A <- generateA("SBM", params)
    
    # E-value with gamma = 0.1
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "ER", "SBM", 2, 2, 0.1, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "ER", "SBM", 2, 2, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej
    
    # E-value with gamma = 0.9
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "ER", "SBM", 2, 2, 0.9, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses SBM with 2 community
    df$time[cnt+3] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=2, s=5, o=100, R=1, loss="l2")$l2.model=="SBM-2")
    )[3]
    df$rej[cnt+3] <- rej
    
    cnt = cnt + length(methods)
    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_K2_delta.RData")
    
  }
  print(delta)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_K2_delta.RData")

df_plot <- df %>% group_by(Method, delta) %>% summarize(rej = mean(rej), time=mean(time))

p1 <- ggplot(df_plot, aes(x=delta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Rejection Rate")+
  theme_bw()
p1

ggsave("~/Documents/Research/network_model_selection/Figures/sbm_K2_delta.pdf", height=4, width=6, unit="in")


p2 <- ggplot(df_plot, aes(x=delta, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Time (sec)")+
  theme_bw()


ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")


## Setting 1 (c) - H0: SBM w/ K=4 vs. H1: SBM w/ K=5
## Increasing community size

niter = 10
n = 1000
alpha = 0.1
beta = 0.8
delta.seq = seq(0.01, 0.20, length=10)
K = 5

B = alpha * (beta*diag(K) + (1-beta)*matrix(1, K, K))

methods = c("e-value0.4", "e-value0.5", "e-value0.6", "NETCROP")

df <- tibble(iter = rep(rep(1:niter,each=length(methods)), length(delta.seq)), 
             Method=rep(methods, niter*length(delta.seq)), 
             delta=0, rej=0, time=0)

cnt = 1


for(delta in delta.seq){
  
  pi = c(rep(1/(K-1) - delta/(K-1), K-1), delta)
  C  = sample(1:K, n, TRUE, pi)
  
  for(iter in 1:niter){
    
    df$delta[cnt:(cnt+length(methods)-1)] <- delta
    
    params = list(B=B, C=C)
    A <- generateA("SBM", params)
    
    # E-value with gamma = 0.1
    # Reject if e-value is greater than 20.
    df$time[cnt] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "SBM", K-1, K, 0.4, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt]  <- rej
    
    # E-value with gamma = 0.5
    # Reject if e-value is greater than 20.
    df$time[cnt+1] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "SBM", K-1, K, 0.5, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+1]  <- rej
    
    # E-value with gamma = 0.9
    # Reject if e-value is greater than 20.
    df$time[cnt+2] <- system.time(rej <- as.numeric(eval_mc(A, "SBM", "SBM", K-1, K, 0.6, nreps = nreps, ncores = detectCores()-1) > 20))[3]
    df$rej[cnt+2]  <- rej
    
    # Reject if it chooses SBM with 5 communities
    df$time[cnt+3] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=2, s=5, o=100, R=1, loss="l2")$l2.model=="SBM-5")
    )[3]
    df$rej[cnt+3] <- rej
    
    cnt = cnt + length(methods)
    save(df, file="~/Documents/Research/network_model_selection/Results/df_sbm_K5_delta.RData")
    
  }
  print(delta)
}

load("~/Documents/Research/network_model_selection/Results/df_sbm_K5_delta.RData")

df_plot <- df %>% group_by(Method, delta) %>% summarize(rej = mean(rej), time=mean(time))

p1 <- ggplot(df_plot, aes(x=delta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Rejection Rate")+
  theme_bw()
p1
ggsave("~/Documents/Research/network_model_selection/Figures/sbm_K5_delta.pdf", height=4, width=6, unit="in")


p2 <- ggplot(df_plot, aes(x=delta, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Time (sec)")+
  theme_bw()

ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")


