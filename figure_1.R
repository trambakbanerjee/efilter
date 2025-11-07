
library(Rfast)
library(foreach)
library(doParallel)
library(adaFilter)
library(ParFilter)

source('funcs.R')

exp4<- function(m,agents,r,all_zero_frac,rho1,rho2,alph,reps){
  
  Sigma_chol_1 = AR1_chol(m,rho = rho1)
  #set.seed(2023)
  Sigma_chol_2 = chol(rho2*matrix(1,agents,agents)+(1-rho2)*diag(agents))
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  result<-foreach(rr = 1:reps,.packages=c("adaFilter","ParFilter"))%dopar%{
    
    source('funcs.R')
    set.seed(rr)
    data <- GenPMat_corr(M = m, Sigma_chol_1,Sigma_chol_2,
                         n = agents, r = r, alternative.frac = 0.01, 
                         all.zero.frac = all_zero_frac,
                         mu = c(-5,-4,-3,3,4,5))
    theta = 1*data$truth.pc
    # AdaFilter
    adaFilter_result <- adaFilter(data$pvalue.mat, r = r,alpha = alph,type.I.err = "FDR")
    adaFilter_fdp = sum((1-theta)*adaFilter_result$decision)/max(sum(adaFilter_result$decision),1)
    adaFilter_etp = sum(theta*adaFilter_result$decision)/max(sum(theta),1)
    
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
    bh_bpvals = bh.func(bpvals,alph)
    bpvals_fdp = sum((1-theta)*bh_bpvals$de)/max(sum(bh_bpvals$de),1)
    bpvals_etp = sum(theta*bh_bpvals$de)/max(sum(theta),1)
    
    ## Cauchy pvals ###############################################
    cauchy_out = cauchy_PCH(pvals,m,agents,r,alph,"FDR")
    cauchy_de = cauchy_out$decisions
    cauchypvals = cauchy_out$pvals
    
    cauchypvals_fdp = sum((1-theta)*cauchy_de)/max(sum(cauchy_de),1)
    cauchypvals_etp = sum(theta*cauchy_de)/max(sum(theta),1)
    #####################################################################
    
    ## e-Filter Bonferroni
    kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
    e_ada_rejs = matrix(0,length(kappa_vals))
    for(k in 1:length(kappa_vals)){
      
      kappa = kappa_vals[k]
      e1 = kappa*(F_j^{kappa-1}) 
      e2 = kappa*(bpvals^{kappa-1}) 
      
      de = e_filter(e2,e1,m,alph,"FDR")$decisions
      e_ada_rejs[k] = sum(de)
    }
    kappa = kappa_vals[which.max(e_ada_rejs)]
    e1 = kappa*(F_j^{kappa-1}) 
    e2 = kappa*(bpvals^{kappa-1}) 
    de = e_filter(e2,e1,m,alph,"FDR")$decisions
    efilter_bonf_fdp = sum((1-theta)*de)/max(sum(de),1)
    efilter_bonf_etp = sum(theta*de)/max(sum(theta),1)
    
    ## e-Filter Cauchy
    kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
    e_ada_rejs = matrix(0,length(kappa_vals))
    for(k in 1:length(kappa_vals)){
      
      kappa = kappa_vals[k]
      e1 = kappa*(cauchy_out$cauchy_eF^{kappa-1}) 
      e2 = kappa*(cauchypvals^{kappa-1}) 
      
      de = e_filter(e2,e1,m,alph,"FDR")$decisions
      e_ada_rejs[k] = sum(de)
    }
    kappa = kappa_vals[which.max(e_ada_rejs)]
    e1 = kappa*(cauchy_out$cauchy_eF^{kappa-1}) 
    e2 = kappa*(cauchypvals^{kappa-1}) 
    de = e_filter(e2,e1,m,alph,"FDR")$decisions
    efilter_cauchy_fdp = sum((1-theta)*de)/max(sum(de),1)
    efilter_cauchy_etp = sum(theta*de)/max(sum(theta),1)
    
    # Naive
    kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
    e_PCH_rejs = matrix(0,length(kappa_vals))
    for(k in 1:length(kappa_vals)){
      
      kappa = kappa_vals[k]
      evals = kappa*(pvals^{kappa-1}) 
      de = e_PCH(evals,m,agents,r,alph,"FDR")$decisions
      e_PCH_rejs[k] = sum(de)
    }
    kappa = kappa_vals[which.max(e_PCH_rejs)]
    evals = kappa*(pvals^{kappa-1}) 
    de = e_PCH(evals,m,agents,r,alph,"FDR")$decisions
    e_PCH_fdp = sum((1-theta)*de)/max(sum(de),1)
    e_PCH_etp = sum(theta*de)/max(sum(theta),1)
    #####################################################
    X_list <- lapply(1:agents, function(x) rnorm(m))
    Rejections <- ParFilter_FDR(p_mat = pvals, X_list=X_list,
                                u = r, q = alph, K = r,
                                method = "Simes", adaptive = TRUE, 
                                cross_weights = FALSE,
                                lambdas = rep(0.50,r))
    de_parfilter = rep(0,m)
    de_parfilter[Rejections] = 1
    parfilter_fdp = sum((1-theta)*de_parfilter)/max(sum(de_parfilter),1)
    parfilter_etp = sum(theta*de_parfilter)/max(sum(theta),1)
    
    #### Cofilter
    fisher_pvals = matrix(0,m,1)
    for(i in 1:m){
      
      cc = pvals[i,]
      cc = cc[order(cc)]
      temp = -2*sum(log(cc[r:agents]))
      fisher_pvals[i] = min(1,1-pchisq(temp,2*(agents-r+1)))
    }
    tau = seq(0.01,1,length.out=50)
    cof_rejs = matrix(0,length(tau),1)
    for(taus in 1:length(tau)){
    
      idx = (fisher_pvals<tau[taus])
      np = (1/tau[taus])*fisher_pvals[idx]
      cof_rejs[taus] = sum(bh.func(np,alph)$de)
    }
    tau_star = tau[which.max(cof_rejs)]
    idx = (fisher_pvals<tau_star)
    np = (1/tau_star)*fisher_pvals[idx]
    cof_rejs = bh.func(np,alph)$de
    de = matrix(0,m,1)
    de[idx] = cof_rejs
    cofilt_fdp = sum((1-theta)*de)/max(sum(de),1)
    cofilt_etp = sum(theta*de)/max(sum(theta),1)
    
    
    return(list('ada_fdp' = adaFilter_fdp,
                'cauchy_fdp' = cauchypvals_fdp,
                'bpvals_fdp' = bpvals_fdp,
                'efilter_bonf_fdp' = efilter_bonf_fdp,
                'efilter_cauchy_fdp' = efilter_cauchy_fdp,
                'e_PCH_fdp' = e_PCH_fdp,
                'parfilter_fdp' = parfilter_fdp,
                'cofilter_fdp' = cofilt_fdp,
                'ada_etp' = adaFilter_etp,
                'cauchy_etp' = cauchypvals_etp,
                'efilter_bonf_etp' = efilter_bonf_etp,
                'efilter_cauchy_etp' = efilter_cauchy_etp,
                'bpvals_etp' = bpvals_etp,
                'e_PCH_etp' =e_PCH_etp,
                'parfilter_etp' = parfilter_etp,
                'cofilter_etp' = cofilt_etp))
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  return(result)
  
}

#-------------------------------------------------------------
m = 10000
alph = 0.2
reps = 500

out8_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.8,alph,reps)
out7_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.7,alph,reps)
out6_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.6,alph,reps)
out5_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.5,alph,reps)
out4_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.4,alph,reps)
out3_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.3,alph,reps)
out2_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.2,alph,reps)
out1_p98 = exp4(m,5,2,all_zero_frac = 0.98,0,0.1,alph,reps)

save.image('figure_1.RData')


ada_fdp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out2_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out3_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out4_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out5_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out6_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out7_p98[[i]]$ada_fdp)),
            mean(sapply(1:reps,function(i) out8_p98[[i]]$ada_fdp)))

ada_etp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out2_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out3_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out4_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out5_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out6_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out7_p98[[i]]$ada_etp)),
            mean(sapply(1:reps,function(i) out8_p98[[i]]$ada_etp)))

parf_fdp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out2_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out3_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out4_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out5_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out6_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out7_p98[[i]]$parfilter_fdp)),
            mean(sapply(1:reps,function(i) out8_p98[[i]]$parfilter_fdp)))
parf_etp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out2_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out3_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out4_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out5_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out6_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out7_p98[[i]]$parfilter_etp)),
             mean(sapply(1:reps,function(i) out8_p98[[i]]$parfilter_etp)))

cof_fdp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out2_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out3_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out4_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out5_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out6_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out7_p98[[i]]$cofilter_fdp)),
             mean(sapply(1:reps,function(i) out8_p98[[i]]$cofilter_fdp)))
cof_etp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out2_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out3_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out4_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out5_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out6_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out7_p98[[i]]$cofilter_etp)),
             mean(sapply(1:reps,function(i) out8_p98[[i]]$cofilter_etp)))

bpvals_fdp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out2_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out3_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out4_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out5_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out6_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out7_p98[[i]]$bpvals_fdp)),
               mean(sapply(1:reps,function(i) out8_p98[[i]]$bpvals_fdp)))

bpvals_etp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out2_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out3_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out4_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out5_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out6_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out7_p98[[i]]$bpvals_etp)),
               mean(sapply(1:reps,function(i) out8_p98[[i]]$bpvals_etp)))

eval_3_fdp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out2_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out3_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out4_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out5_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out6_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out7_p98[[i]]$efilter_bonf_fdp)),
               mean(sapply(1:reps,function(i) out8_p98[[i]]$efilter_bonf_fdp)))

eval_3_etp = c(mean(sapply(1:reps,function(i) out1_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out2_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out3_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out4_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out5_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out6_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out7_p98[[i]]$efilter_bonf_etp)),
               mean(sapply(1:reps,function(i) out8_p98[[i]]$efilter_bonf_etp)))


library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1 = data.frame('fdp'=c(ada_fdp,parf_fdp,cof_fdp,bpvals_fdp,eval_3_fdp),
                       'etp' = c(ada_etp,parf_etp,cof_etp,bpvals_etp,eval_3_etp),
                       'rho' = rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),5),
                       'method' = rep(c('AdaFilter','ParFilter','CoFilter',
                                        'Bonferroni','e-Filter'),
                                      each = 8))

g1<-ggplot(data=plotdata1,aes(x=rho,y=fdp,col=factor(method),shape=factor(method)))+
  geom_line()+geom_point()+geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  ylab('FDR')+theme_bw()+xlab(TeX('$\\rho$'))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

g2<-ggplot(data=plotdata1,aes(x=rho,y=etp,col=factor(method),shape=factor(method)))+
  geom_line()+geom_point()+
  ylab('Recall')+theme_bw()+xlab(TeX('$\\rho$'))+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggarrange(g1,g2,ncol=2,nrow=1,common.legend = TRUE)


