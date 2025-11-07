library(adaFilter)
library(Rfast)
library(ggplot2)
library(ggpubr)
library(latex2exp)

source('funcs.R')

m = 10000
alph = 1
all_zero_frac = 0.98
rho1 = 0
rho2 = 0.7
agents = 2
r = 2
rr = 1

Sigma_chol_1 = AR1_chol(m,rho = rho1)
#set.seed(2023)
Sigma_chol_2 = chol(rho2*matrix(1,agents,agents)+(1-rho2)*diag(agents))


set.seed(rr)
data <- GenPMat_corr(M = m, Sigma_chol_1,Sigma_chol_2,
                     n = agents, r = r, alternative.frac = 0.01, 
                     all.zero.frac = all_zero_frac,
                     mu = c(-6,-5,-4,4,5,6))
theta = 1*data$truth.pc
# AdaFilter
adaFilter_result <- adaFilter(data$pvalue.mat, r = r,alpha = alph,type.I.err = "PFER")
adaFilter_fdp = sum((1-theta)*adaFilter_result$decision)
adaFilter_etp = sum(theta*adaFilter_result$decision)/max(sum(theta),1)
adaFilter_gam_1 = max(adaFilter_result$selection.p[adaFilter_result$decision==1])

pvals = data$pvalue.mat

## Bonferroni pvals
bpvals = matrix(0,m,1)
F_j = matrix(0,m,1)
for(i in 1:m){
  
  cc = pvals[i,]
  cc = cc[order(cc)]
  bpvals[i] = min(1,(agents-r+1)*cc[r])
  F_j[i] = min(1,(agents-r+1)*cc[r-1])
}
bh_bpvals = 1*(bpvals<=(alph/m))
bpvals_fdp = sum((1-theta)*bh_bpvals)
bpvals_etp = sum(theta*bh_bpvals)/max(sum(theta),1)

## E-Filter
kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
e_ada_rejs = matrix(0,length(kappa_vals))
for(k in 1:length(kappa_vals)){
  
  kappa = kappa_vals[k]
  e1 = kappa*(F_j^{kappa-1}) 
  e2 = kappa*(bpvals^{kappa-1}) 
  
  ff<-function(gam,e1,e2,alpha){
    
    return(gam*sum(1*(e1>(1/gam)))-alpha)
  }
  
  sol = try(uniroot(ff,c(0,alph),e1=e1,e2=e2,alpha=alph))
  gam_BH = if(inherits(sol,"try-error")) NA else sol$root
  
  if(!is.na(gam_BH)){
    de = 1*(e2>(1/gam_BH))
    e_ada_rejs[k] = sum(de)
  }
  # print(gam_BH)
}
kappa = kappa_vals[which.max(e_ada_rejs)]
e1 = kappa*(F_j^{kappa-1}) 
e2 = kappa*(bpvals^{kappa-1}) 

ff<-function(gam,e1,e2,alpha){
  
  return(gam*sum(1*(e1>(1/gam)))-alpha)
}
gam_BH_1 = uniroot(ff,c(0,alph),e1=e1,e2=e2,alpha=alph)$root

de = 1*(e2>(1/gam_BH_1))
e_fdp_3 = sum((1-theta)*de)
e_etp_3 = sum(theta*de)/max(sum(theta),1)

adaFilter_gam_1
gam_BH_1


p1 = -log(pvals[theta==0,1],base=10)
p2 = -log(pvals[theta==0,2],base=10)
p_ada = -log(pvals[theta==0 & adaFilter_result$decision>0,],base=10)
p_eada = -log(pvals[theta==0 & de>0,],base=10)

plotdata1 = data.frame('p1'=p1,'p2'=p2,
                       'p1_ada' =c(p_ada[,1],rep(NA,length(p1)-nrow(p_ada))),
                       'p2_ada'=c(p_ada[,2],rep(NA,length(p1)-nrow(p_ada))))

g11<-ggplot()+geom_point(data=plotdata1,aes(x=p1,y=p2),col='black',size = 1)+
  geom_point(data=plotdata1,aes(x=p1_ada,y=p2_ada),shape=22,
             col='red',fill='red',size = 2)+
  geom_rect(aes(xmin = -log((kappa*gam_BH_1)^{1/(1-kappa)},base=10), xmax = Inf, ymin = -log((kappa*gam_BH_1)^{1/(1-kappa)},base=10),ymax=Inf), 
            color = "black",fill=NA,linetype=2)+
  geom_rect(aes(xmin = -log(adaFilter_gam_1,base=10), xmax=Inf,ymin=-log(adaFilter_gam_1,base=10),ymax=Inf), 
            color = "blue",fill=NA)+
  #ylim(2,5)+xlim(2,5)+
  annotate("text", x = 8, y = 10, label = TeX('$\\rho=0.7$')) +
  ylab(TeX('$-\\log_{10}(p_{2j})$ when $H_{0j}^{2/2}$ is true'))+xlab(TeX('$-\\log_{10}(p_{1j})$ when $H_{0j}^{2/2}$ is true'))+
  theme_bw()

p_bonf_ef = -log(pvals[theta==1 & bh_bpvals>0 & de>0,],base=10)
p_ef_only = -log(pvals[theta==1 & de>0 & bh_bpvals==0,],base=10)

plotdata1 = data.frame('p1_bonf' =c(p_bonf_ef[,1]),
                       'p2_bonf'=c(p_bonf_ef[,2]),
                       'p1_ef' =c(p_ef_only[,1],rep(NA,nrow(p_bonf_ef)-nrow(p_ef_only))),
                       'p2_ef'=c(p_ef_only[,2],rep(NA,nrow(p_bonf_ef)-nrow(p_ef_only))))

g21<- ggplot()+geom_point(data=plotdata1,aes(x=p1_bonf,y=p2_bonf),color='black',size = 1)+
  geom_point(data=plotdata1,aes(x=p1_ef,y=p2_ef),shape=22,
             col='red',fill='red',size = 2)+
  geom_rect(aes(xmin = -log(alph/m,base=10), xmax=Inf,
                ymin=-log(alph/m,base=10),ymax=Inf), 
            color = "blue",fill=NA)+
  geom_rect(aes(xmin = -log((kappa*gam_BH_1)^{1/(1-kappa)},base=10), xmax = Inf, 
                ymin = -log((kappa*gam_BH_1)^{1/(1-kappa)},base=10),ymax=Inf), 
            color = "black",fill=NA,linetype=2)+
  annotate("text", x = 10, y = 13, label = TeX('$\\rho=0.7$')) +
  ylab(TeX('$-\\log_{10}(p_{2j})$ when $H_{0j}^{2/2}$ is false'))+xlab(TeX('$-\\log_{10}(p_{1j})$ when $H_{0j}^{2/2}$ is false'))+
  theme_bw()
  
  
  #--------------------------------------------------------------
rho2 = 0.9
Sigma_chol_2 = chol(rho2*matrix(1,agents,agents)+(1-rho2)*diag(agents))


set.seed(rr)
data <- GenPMat_corr(M = m, Sigma_chol_1,Sigma_chol_2,
                     n = agents, r = r, alternative.frac = 0.01, 
                     all.zero.frac = all_zero_frac,
                     mu = c(-6,-5,-4,4,5,6))
theta = 1*data$truth.pc
# AdaFilter
adaFilter_result <- adaFilter(data$pvalue.mat, r = r,alpha = alph,type.I.err = "PFER")
adaFilter_fdp = sum((1-theta)*adaFilter_result$decision)
adaFilter_etp = sum(theta*adaFilter_result$decision)/max(sum(theta),1)
adaFilter_gam_2 = max(adaFilter_result$selection.p[adaFilter_result$decision==1])

pvals = data$pvalue.mat

## Bonferroni pvals
bpvals = matrix(0,m,1)
F_j = matrix(0,m,1)
for(i in 1:m){
  
  cc = pvals[i,]
  cc = cc[order(cc)]
  bpvals[i] = min(1,(agents-r+1)*cc[r])
  F_j[i] = min(1,(agents-r+1)*cc[r-1])
}
bh_bpvals = 1*(bpvals<=(alph/m))
bpvals_fdp = sum((1-theta)*bh_bpvals)
bpvals_etp = sum(theta*bh_bpvals)/max(sum(theta),1)

## E-Filter
kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
e_ada_rejs = matrix(0,length(kappa_vals))
for(k in 1:length(kappa_vals)){
  
  kappa = kappa_vals[k]
  e1 = kappa*(F_j^{kappa-1}) 
  e2 = kappa*(bpvals^{kappa-1}) 
  
  ff<-function(gam,e1,e2,alpha){
    
    return(gam*sum(1*(e1>(1/gam)))-alpha)
  }
  
  sol = try(uniroot(ff,c(0,alph),e1=e1,e2=e2,alpha=alph))
  gam_BH = if(inherits(sol,"try-error")) NA else sol$root

  if(!is.na(gam_BH)){
    de = 1*(e2>(1/gam_BH))
    e_ada_rejs[k] = sum(de)
  }
  # print(gam_BH)
}
kappa = kappa_vals[which.max(e_ada_rejs)]
e1 = kappa*(F_j^{kappa-1}) 
e2 = kappa*(bpvals^{kappa-1}) 

ff<-function(gam,e1,e2,alpha){
  
  return(gam*sum(1*(e1>(1/gam)))-alpha)
}
gam_BH_2 = uniroot(ff,c(0,alph),e1=e1,e2=e2,alpha=alph)$root

de = 1*(e2>(1/gam_BH_2))
e_fdp_3 = sum((1-theta)*de)
e_etp_3 = sum(theta*de)/max(sum(theta),1)

adaFilter_gam_2
gam_BH_2


#-------------------------------------------------------------

p1 = -log(pvals[theta==0,1],base=10)
p2 = -log(pvals[theta==0,2],base=10)
p_ada = -log(pvals[theta==0 & adaFilter_result$decision>0,],base=10)
p_eada = -log(pvals[theta==0 & de>0,],base=10)

plotdata2 = data.frame('p1'=p1,'p2'=p2,
                       'p1_ada' =c(p_ada[,1],rep(NA,length(p1)-nrow(p_ada))),
                       'p2_ada'=c(p_ada[,2],rep(NA,length(p1)-nrow(p_ada))),
                       'p1_eada' =c(p_eada[1],rep(NA,length(p1)-1)),
                       'p2_eada'=c(p_eada[2],rep(NA,length(p1)-1)))

g12<-ggplot()+geom_point(data=plotdata2,aes(x=p1,y=p2),col='black',size = 1)+
  geom_point(data=plotdata2,aes(x=p1_ada,y=p2_ada),shape=22,
             col='red',fill='red',size = 2)+
  geom_rect(aes(xmin = -log((kappa*gam_BH_2)^{1/(1-kappa)},base=10), xmax = Inf, ymin = -log((kappa*gam_BH_2)^{1/(1-kappa)},base=10),ymax=Inf), 
            color = "black",fill=NA,linetype=2)+
  geom_rect(aes(xmin = -log(adaFilter_gam_2,base=10), xmax=Inf,ymin=-log(adaFilter_gam_2,base=10),ymax=Inf), 
            color = "blue",fill=NA)+
  annotate("text", x = 8, y = 10, label = TeX('$\\rho=0.9$')) +
  ylab(TeX('$-\\log_{10}(p_{2j})$ when $H_{0j}^{2/2}$ is true'))+xlab(TeX('$-\\log_{10}(p_{1j})$ when $H_{0j}^{2/2}$ is true'))+
  theme_bw()

p_bonf_ef = -log(pvals[theta==1 & bh_bpvals>0 & de>0,],base=10)
p_ef_only = -log(pvals[theta==1 & de>0 & bh_bpvals==0,],base=10)

plotdata2 = data.frame('p1_bonf' =c(p_bonf_ef[,1]),
                       'p2_bonf'=c(p_bonf_ef[,2]),
                       'p1_ef' =c(p_ef_only[,1],rep(NA,nrow(p_bonf_ef)-nrow(p_ef_only))),
                       'p2_ef'=c(p_ef_only[,2],rep(NA,nrow(p_bonf_ef)-nrow(p_ef_only))))

g22<- ggplot()+geom_point(data=plotdata2,aes(x=p1_bonf,y=p2_bonf),color='black',size = 1)+
  geom_point(data=plotdata2,aes(x=p1_ef,y=p2_ef),shape=22,
             col='red',fill='red',size = 2)+
  geom_rect(aes(xmin = -log(alph/m,base=10), xmax=Inf,
                ymin=-log(alph/m,base=10),ymax=Inf), 
            color = "blue",fill=NA)+
  geom_rect(aes(xmin = -log((kappa*gam_BH_2)^{1/(1-kappa)},base=10), xmax = Inf, 
                ymin = -log((kappa*gam_BH_2)^{1/(1-kappa)},base=10),ymax=Inf), 
            color = "black",fill=NA,linetype=2)+
  annotate("text", x = 10, y = 13, label = TeX('$\\rho=0.9$')) +
  ylab(TeX('$-\\log_{10}(p_{2j})$ when $H_{0j}^{2/2}$ is false'))+
  xlab(TeX('$-\\log_{10}(p_{1j})$ when $H_{0j}^{2/2}$ is false'))+
  theme_bw()

ggarrange(g11,g12,ncol=2)

ggarrange(g21,g22,ncol=2)


