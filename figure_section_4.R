
library(MASS)
library(mvtnorm)
source('funcs.R')


rho.vals = c(-0.3,0,0.3,0.5,0.7,0.9,0.99)

set.seed(1)
analyze<-function(mu,rho.vals){
  
  x = seq(0,1,length.out=50)
  x = x[-c(1,50)]
  means = mu
  
  d1 = d2 = k = matrix(0,length(rho.vals),1)
  
  for(j in 1:length(rho.vals)){
    
    rho =rho.vals[j]
    
    Sigma = rbind(c(1,rho,rho,rho),
                  c(rho,1,rho,rho),
                  c(rho,rho,1,rho),
                  c(rho,rho,rho,1)) #
    tau = x
    
    C2 = matrix(0,length(x),1)
    for(i in 1:length(x)){
      ztau = qnorm(1-0.5*tau[i])
      C2[i] = 1-(pmvnorm(lower=c(-ztau,-ztau,-ztau,-ztau),
                         upper=c(ztau,ztau,ztau,ztau),
                         mean = means,
                         sigma = Sigma))
    }
    yy2 = function(Cval,d){
      x = seq(0,1,length.out=50)
      x = x[-c(1,50)]
      return(min(Cval-x^{d}))
    }
    d2[j] = uniroot(yy2, c(0, 1), tol = 0.0001, Cval = C2)$root
    
    C1 = matrix(0,length(x),1)
    for(i in 1:length(x)){
      ztau = qnorm(1-0.5*tau[i])
      C1[i] = pmvnorm(lower=c(ztau,ztau,ztau,ztau),
                      upper=c(Inf,Inf,Inf,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,ztau,-Inf,ztau),
                upper=c(Inf,Inf,-ztau,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,-Inf,ztau,ztau),
                upper=c(Inf,-ztau,Inf,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,-Inf,-Inf,ztau),
                upper=c(Inf,-ztau,-ztau,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,-Inf,-Inf,ztau),
                upper=c(-ztau,-ztau,-ztau,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,-Inf,ztau,ztau),
                upper=c(-ztau,-ztau,Inf,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,ztau,ztau,ztau),
                upper=c(-ztau,Inf,Inf,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,ztau,-Inf,ztau),
                upper=c(-ztau,Inf,-ztau,Inf),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,ztau,ztau,-Inf),
                upper=c(Inf,Inf,Inf,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,ztau,-Inf,-Inf),
                upper=c(Inf,Inf,-ztau,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,-Inf,ztau,-Inf),
                upper=c(Inf,-ztau,Inf,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(ztau,-Inf,-Inf,-Inf),
                upper=c(Inf,-ztau,-ztau,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,-Inf,-Inf,-Inf),
                upper=c(-ztau,-ztau,-ztau,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,-Inf,ztau,-Inf),
                upper=c(-ztau,-ztau,Inf,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,ztau,ztau,-Inf),
                upper=c(-ztau,Inf,Inf,ztau),mean = means,sigma = Sigma)+
        pmvnorm(lower=c(-Inf,ztau,-Inf,-Inf),
                upper=c(-ztau,Inf,-ztau,ztau),mean = means,sigma = Sigma)
    }
    yy1 = function(Cval,d){
      x = seq(0,1,length.out=50)
      x = x[-c(1,50)]
      return(min(x^{d}-Cval))
    }
    d1[j] = uniroot(yy1, c(1, 10), tol = 0.0001, Cval = C1)$root
    
    k[j] = max(0,d2[j]+1-d1[j])
  }
  
  return(data.frame('rho'=rho.vals,'y'=d1-d2+k,'k'=k))
}

least_fav = analyze(c(0,3,3,3),rho.vals)
one_nonnull = analyze(c(0,3,0,0),rho.vals)
two_nonnull = analyze(c(0,3,0,3),rho.vals)
global_null = analyze(c(0,0,0,0),rho.vals)

library(ggplot2)
library(ggpubr)
library(latex2exp)


g4<- ggplot(data=global_null)+
  geom_point(aes(x=rho,y=k),col='red',shape=15)+
  geom_line(aes(x=rho,y=k))+
  geom_hline(yintercept = 1,linetype='dashed')+
  scale_x_continuous(breaks=rho.vals)+
  ggtitle('Global Null') +
  ylab(TeX('$\\kappa^{*}$'))+xlab(TeX('$\\rho$'))+
  theme_bw()

g2<- ggplot(data=one_nonnull)+
  geom_point(aes(x=rho,y=k),col='red',shape=15)+
  ylim(c(0,0.05))+
  geom_line(aes(x=rho,y=k))+
  #geom_hline(yintercept = 1,linetype='dashed')+
  scale_x_continuous(breaks=rho.vals)+
  ggtitle('One Non-Null') +
  ylab(TeX('$\\kappa^{*}$'))+xlab(TeX('$\\rho$'))+
  theme_bw()

g3<- ggplot(data=two_nonnull)+
  geom_point(aes(x=rho,y=k),col='red',shape=15)+
  geom_line(aes(x=rho,y=k))+
  ylim(c(0,0.05))+
  #geom_hline(yintercept = 1,linetype='dashed')+
  scale_x_continuous(breaks=rho.vals)+
  ggtitle('Two Non-Nulls') +
  ylab(TeX('$\\kappa^{*}$'))+xlab(TeX('$\\rho$'))+
  theme_bw()

g1<-ggplot(data=least_fav)+
  geom_point(aes(x=rho,y=k),col='red',shape=15)+
  geom_line(aes(x=rho,y=k))+
  ylim(c(0,0.05))+
  #geom_hline(yintercept = 1,linetype='dashed')+
  scale_x_continuous(breaks=rho.vals)+
  ggtitle('Least Favorable Null') +
  ylab(TeX('$\\kappa^{*}$'))+xlab(TeX('$\\rho$'))+
  theme_bw()

ggarrange(g2,g3,g1,g4,ncol=2,nrow=2)


