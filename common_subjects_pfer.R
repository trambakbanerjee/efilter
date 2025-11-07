
#Common case

library(adaFilter)
library(Rfast)
library(foreach)
library(doParallel)
library(adaFilter)

source('funcs.R')

exp4<- function(m,r,q,all_zero_frac,alph,reps){
  
  agents = 5
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  result<-foreach(rr = 1:reps,.packages=c("adaFilter"))%dopar%{
    
    source('funcs.R')
    set.seed(rr)
    data <- GenPMat_commoncases_n5(M = m, r = r,s = 5000,
                                   q = q,
                         alternative.frac = 0.01, 
                         all.zero.frac = all_zero_frac,
                         mu = c(-4,-3,3,4),
                         one.sided = FALSE)
    theta = 1*data$truth.pc
    # AdaFilter
    adaFilter_result <- adaFilter(data$pvalue.mat, r = r,alpha = alph,type.I.err = "PFER")
    adaFilter_fdp = sum((1-theta)*adaFilter_result$decision)
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
    bh_bpvals = 1*(bpvals<=(alph/m))
    bpvals_fdp = sum((1-theta)*bh_bpvals)
    bpvals_etp = sum(theta*bh_bpvals)/max(sum(theta),1)
    
    ## Cauchy pvals ###############################################
    cauchy_out = cauchy_PCH(pvals,m,agents,r,alph,"PFER")
    cauchy_de = cauchy_out$decisions
    cauchypvals = cauchy_out$pvals
    
    cauchypvals_fdp = sum((1-theta)*cauchy_de)
    cauchypvals_etp = sum(theta*cauchy_de)/max(sum(theta),1)
    #####################################################################
    
    ## e-Filter Bonferroni
    kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
    e_ada_rejs = matrix(0,length(kappa_vals))
    for(k in 1:length(kappa_vals)){
      
      kappa = kappa_vals[k]
      e1 = kappa*(F_j^{kappa-1}) 
      e2 = kappa*(bpvals^{kappa-1}) 
      
      de = e_filter(e2,e1,m,alph,"PFER")$decisions
      e_ada_rejs[k] = sum(de)
    }
    kappa = kappa_vals[which.max(e_ada_rejs)]
    e1 = kappa*(F_j^{kappa-1}) 
    e2 = kappa*(bpvals^{kappa-1}) 
    de = e_filter(e2,e1,m,alph,"PFER")$decisions
    efilter_bonf_fdp = sum((1-theta)*de)
    efilter_bonf_etp = sum(theta*de)/max(sum(theta),1)
    
    ## e-Filter Cauchy
    kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
    e_ada_rejs = matrix(0,length(kappa_vals))
    for(k in 1:length(kappa_vals)){
      
      kappa = kappa_vals[k]
      e1 = kappa*(cauchy_out$cauchy_eF^{kappa-1}) 
      e2 = kappa*(cauchypvals^{kappa-1}) 
      
      de = e_filter(e2,e1,m,alph,"PFER")$decisions
      e_ada_rejs[k] = sum(de)
    }
    kappa = kappa_vals[which.max(e_ada_rejs)]
    e1 = kappa*(cauchy_out$cauchy_eF^{kappa-1}) 
    e2 = kappa*(cauchypvals^{kappa-1}) 
    de = e_filter(e2,e1,m,alph,"PFER")$decisions
    efilter_cauchy_fdp = sum((1-theta)*de)
    efilter_cauchy_etp = sum(theta*de)/max(sum(theta),1)
    
    # Naive
    kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
    e_PCH_rejs = matrix(0,length(kappa_vals))
    for(k in 1:length(kappa_vals)){
      
      kappa = kappa_vals[k]
      evals = kappa*(pvals^{kappa-1}) 
      de = e_PCH(evals,m,agents,r,alph,"PFER")$decisions
      e_PCH_rejs[k] = sum(de)
    }
    kappa = kappa_vals[which.max(e_PCH_rejs)]
    evals = kappa*(pvals^{kappa-1}) 
    de = e_PCH(evals,m,agents,r,alph,"PFER")$decisions
    e_PCH_fdp = sum((1-theta)*de)
    e_PCH_etp = sum(theta*de)/max(sum(theta),1)
    #####################################################
 
    return(list('ada_fdp' = adaFilter_fdp,
                'cauchy_fdp' = cauchypvals_fdp,
                'bpvals_fdp' = bpvals_fdp,
                'efilter_bonf_fdp' = efilter_bonf_fdp,
                'efilter_cauchy_fdp' = efilter_cauchy_fdp,
                'e_PCH_fdp' = e_PCH_fdp,
                'ada_etp' = adaFilter_etp,
                'cauchy_etp' = cauchypvals_etp,
                'efilter_bonf_etp' = efilter_bonf_etp,
                'efilter_cauchy_etp' = efilter_cauchy_etp,
                'bpvals_etp' = bpvals_etp,
                'e_PCH_etp' =e_PCH_etp))
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  return(result)
  
}

#-------------------------------------------------------------
m = 10000
alph = 1
reps = 500

r = 2
outa_r2 = exp4(m,r,100,all_zero_frac = 0.98,alph,reps)
outb_r2 = exp4(m,r,300,all_zero_frac = 0.98,alph,reps)
outc_r2 = exp4(m,r,500,all_zero_frac = 0.98,alph,reps)
outd_r2 = exp4(m,r,700,all_zero_frac = 0.98,alph,reps)
oute_r2 = exp4(m,r,900,all_zero_frac = 0.98,alph,reps)

save.image('common_subjects_pfer.RData')


########################################################
result = list('list1'=outa_r2,'list2'=outb_r2,'list3'=outc_r2,
              'list4'=outd_r2,'list5'=oute_r2)

ada_fdp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
ada_fdp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
ada_etp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
ada_etp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))

bonf_fdp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
bonf_fdp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
bonf_etp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
bonf_etp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))

efilter_bonf_fdp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
efilter_bonf_fdp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
efilter_bonf_etp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
efilter_bonf_etp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))

efilter_cauchy_fdp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
efilter_cauchy_fdp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
efilter_cauchy_etp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
efilter_cauchy_etp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))

cauchy_fdp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
cauchy_fdp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
cauchy_etp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
cauchy_etp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))

epch_fdp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
epch_fdp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
epch_etp = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
epch_etp_sd = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))

library(ggplot2)
library(ggpubr)
library(latex2exp)

plotdata1 = data.frame('fdp'=c(ada_fdp,bonf_fdp,cauchy_fdp,
                               efilter_bonf_fdp,
                               efilter_cauchy_fdp,
                               epch_fdp),
                       'fdp_sd'=c(ada_fdp_sd,
                                  bonf_fdp_sd,cauchy_fdp_sd,
                                  efilter_bonf_fdp_sd,
                                  efilter_cauchy_fdp_sd,
                                  epch_fdp_sd),
                       'etp' = c(ada_etp,bonf_etp,cauchy_etp,
                                 efilter_bonf_etp,
                                 efilter_cauchy_etp,
                                 epch_etp),
                       'etp_sd' = c(ada_etp_sd,
                                    bonf_etp_sd,cauchy_etp_sd,
                                    efilter_bonf_etp_sd,
                                    efilter_cauchy_etp_sd,
                                    epch_etp_sd),
                       'q' = rep(c(100,300,500,700,900),6),
                       'method' = rep(c('AdaFilter', 'Bonferroni','Cauchy','e-Filter B',
                                        'e-Filter C','e-PCH'),
                                      each = 5))
library(ggplot2)
library(ggpubr)
g1<-ggplot(data=plotdata1,aes(x=q,y=fdp,col=factor(method),shape=factor(method)))+
  geom_line(size=1)+geom_point(size=2)+
  geom_hline(yintercept = alph,linetype='dotted', size=1,col = 'black')+
  geom_errorbar(aes(ymin = fdp - fdp_sd, 
                    ymax = fdp + fdp_sd),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  ylab('PFER')+theme_bw()+xlab('q')+
  scale_x_continuous(breaks=c(100,300,500,700,900))+
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

g2<-ggplot(data=plotdata1,aes(x=q,y=etp,col=factor(method),shape=factor(method)))+
  geom_line(size=1)+geom_point(size=2)+
  geom_errorbar(aes(ymin = etp - etp_sd, 
                    ymax = etp + etp_sd),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  ylab('Recall')+theme_bw()+xlab('q')+
  scale_x_continuous(breaks=c(100,300,500,700,900))+
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

