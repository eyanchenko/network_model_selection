## Data fission and UI for model selection on networks
## Eric Yanchenko
## May 1, 2025

source("~/Documents/Research/network_model_selection/netcrop.R")
source("~/Documents/Research/network_model_selection/functions.R")

## Setting 1 - H0: SBM w/ K=1 (ER) vs. H1: SBM w/ K=2
niter = 10
n = 1000
a = 0.05
delta.seq = seq(0, 0.1, length=11)
pp = 0.5

df <- tibble(iter = rep(rep(1:niter,each=2), length(delta.seq)), 
             Method=rep(c("e-value", "NETCROP"), niter*length(delta.seq)), 
             delta=0, rej=0, time=0)

cnt = 1
for(delta in delta.seq){
  B = matrix(c(a+delta,a,a,a+delta), ncol=2)
  prop = 0.75
  C = c(rep(1,prop*n), rep(2,(1-prop)*n))
  
  for(iter in 1:niter){
    
    df$delta[c(cnt,cnt+1)] <- delta
    
    params = list(B=B, C=C)
    A <- generateA("SBM", params)

    
    # Reject if e-value is greater than 10.
    df$time[cnt] <- system.time(rej <- as.numeric(eval(A, "ER", "SBM", 2, 2, pp) > 10))[3]
    df$rej[cnt]  <- rej
    
    # Reject if it chooses SBM with 2 community
    df$time[cnt+1] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=2, s=5, o=100, R=1, loss="l2")$l2.model=="SBM-2")
    )[3]
    df$rej[cnt+1] <- rej
    
    cnt = cnt + 2
  }
  print(delta)
}


df_plot <- df %>% group_by(Method, delta) %>% summarize(rej = mean(rej), time=mean(time))

p1 <- ggplot(df_plot, aes(x=delta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Rejection Rate")+
  theme_bw()

p2 <- ggplot(df_plot, aes(x=delta, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Time (sec)")+
  theme_bw()

library(ggpubr)
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")


# Our method does a lot better.
# We are a little slower but nothing major





## Setting 2 - H0: SBM w/ K=2  vs. H1: DCBM w/ K=2
niter = 10
n = 1000
a = 0.10
b = 0.20
pp = 0.5
delta.seq <- seq(0, 0.5, length=21)

df <- tibble(iter = rep(rep(1:niter,each=2), length(delta.seq)), 
             Method=rep(c("e-value", "NETCROP"), niter*length(delta.seq)), 
             delta=0, rej=0, time=0)

cnt = 1
for(delta in delta.seq){
  B = matrix(c(a+b,a,a,a+b), ncol=2)
  prop = 0.75
  C = c(rep(1,prop*n), rep(2,(1-prop)*n))

  
  for(iter in 1:niter){
    
    df$delta[c(cnt,cnt+1)] <- delta
    
    psi <- runif(n, 0.5 - delta/2, 0.5 + delta/2)
    params = list(psi=psi, B=B, C=C)
    A <- generateA("DCBM", params)
    
    
    # Reject if e-value is greater than 10.
    df$time[cnt] <- system.time(rej <- as.numeric(eval(A, "SBM", "DCBM", 2, 2, pp) > 10))[3]
    df$rej[cnt]  <- rej
    
    # Reject if it chooses DCBM with 2 community
    df$time[cnt+1] <- system.time(
      rej <- as.numeric(croissant.blockmodel(A=A, K.CAND=2, s=5, o=100, R=1, loss="l2")$l2.model=="DCSBM-2")
    )[3]
    df$rej[cnt+1] <- rej
    
    cnt = cnt + 2
  }
  print(delta)
}


df_plot <- df %>% group_by(Method, delta) %>% summarize(rej = mean(rej), time=mean(time))

p1 <- ggplot(df_plot, aes(x=delta, y=rej, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Rejection Rate")+
  theme_bw()

p2 <- ggplot(df_plot, aes(x=delta, y=time, color=Method))+
  geom_point()+
  geom_line()+
  xlab(expression(delta))+
  ylab("Time (sec)")+
  theme_bw()
ggarrange(p1, p2, ncol=2, common.legend = TRUE, legend="bottom")

# Our method does a lot better in this setting as well
# Can also vary n 

# Next steps:
# Make sure how I'm fitting the DCBM is kosher (just with spectral clustering)

# Test H0: DCBM w/ K=2 vs. H1: PABM w/ K=2
# Add Larger number of communities
# Add RDPG, latent space dimensions


