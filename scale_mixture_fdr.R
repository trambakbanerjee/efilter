
#Scale mixture of Gaussians

library(adaFilter)
library(Rfast)
library(foreach)
library(doParallel)
library(adaFilter)

source('funcs.R')

exp4<- function(m,agents,r,all_zero_frac,rho1,rho2,alph,reps){
  
  Sigma_chol_1 = AR1_chol(m,rho = rho1)
  #set.seed(2023)
  Sigma_chol_2 = chol(rho2*matrix(1,agents,agents)+(1-rho2)*diag(agents))
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  result<-foreach(rr = 1:reps,.packages=c("adaFilter"))%dopar%{
    
    source('funcs.R')
    set.seed(rr)
    data <- GenPMat_corr_scalemix(M = m, Sigma_chol_1,Sigma_chol_2,
                         n = agents, r = r, alternative.frac = 0.01, 
                         all.zero.frac = all_zero_frac,
                         mu = c(-6,-5,-4,4,5,6))
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
alph = 0.2
reps = 500

outa_p98 = exp4(m,2,2,all_zero_frac = 0.98,0.5,0.8,alph,reps)
outb_p98 = exp4(m,4,2,all_zero_frac = 0.98,0.5,0.8,alph,reps)
outc_p98 = exp4(m,8,2,all_zero_frac = 0.98,0.5,0.8,alph,reps)
outd_p98 = exp4(m,4,4,all_zero_frac = 0.98,0.5,0.8,alph,reps)
oute_p98 = exp4(m,8,4,all_zero_frac = 0.98,0.5,0.8,alph,reps)
outf_p98 = exp4(m,8,8,all_zero_frac = 0.98,0.5,0.8,alph,reps)

out1_p98 = exp4(m,2,2,all_zero_frac = 0.98,0.5,0.6,alph,reps)
out2_p98 = exp4(m,4,2,all_zero_frac = 0.98,0.5,0.6,alph,reps)
out3_p98 = exp4(m,8,2,all_zero_frac = 0.98,0.5,0.6,alph,reps)
out4_p98 = exp4(m,4,4,all_zero_frac = 0.98,0.5,0.6,alph,reps)
out5_p98 = exp4(m,8,4,all_zero_frac = 0.98,0.5,0.6,alph,reps)
out6_p98 = exp4(m,8,8,all_zero_frac = 0.98,0.5,0.6,alph,reps)

out7_p98 = exp4(m,2,2,all_zero_frac = 0.98,0.5,0.4,alph,reps)
out8_p98 = exp4(m,4,2,all_zero_frac = 0.98,0.5,0.4,alph,reps)
out9_p98 = exp4(m,8,2,all_zero_frac = 0.98,0.5,0.4,alph,reps)
out10_p98 = exp4(m,4,4,all_zero_frac = 0.98,0.5,0.4,alph,reps)
out11_p98 = exp4(m,8,4,all_zero_frac = 0.98,0.5,0.4,alph,reps)
out12_p98 = exp4(m,8,8,all_zero_frac = 0.98,0.5,0.4,alph,reps)

out13_p98 = exp4(m,2,2,all_zero_frac = 0.98,0.5,0.2,alph,reps)
out14_p98 = exp4(m,4,2,all_zero_frac = 0.98,0.5,0.2,alph,reps)
out15_p98 = exp4(m,8,2,all_zero_frac = 0.98,0.5,0.2,alph,reps)
out16_p98 = exp4(m,4,4,all_zero_frac = 0.98,0.5,0.2,alph,reps)
out17_p98 = exp4(m,8,4,all_zero_frac = 0.98,0.5,0.2,alph,reps)
out18_p98 = exp4(m,8,8,all_zero_frac = 0.98,0.5,0.2,alph,reps)

out19_p98 = exp4(m,2,2,all_zero_frac = 0.98,0.5,0,alph,reps)
out20_p98 = exp4(m,4,2,all_zero_frac = 0.98,0.5,0,alph,reps)
out21_p98 = exp4(m,8,2,all_zero_frac = 0.98,0.5,0,alph,reps)
out22_p98 = exp4(m,4,4,all_zero_frac = 0.98,0.5,0,alph,reps)
out23_p98 = exp4(m,8,4,all_zero_frac = 0.98,0.5,0,alph,reps)
out24_p98 = exp4(m,8,8,all_zero_frac = 0.98,0.5,0,alph,reps)

save.image('scale_mixture_fdr.RData')


########################################################
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)

ada_fdp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp))))
ada_etp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp))))

bonf_fdp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp))))
bonf_etp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp))))

efilter_bonf_fdp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp))))
efilter_bonf_etp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp))))

efilter_cauchy_fdp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp))))
efilter_cauchy_etp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp))))

cauchy_fdp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp))))
cauchy_etp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp))))

epch_fdp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp))))
epch_etp = mean(sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp))))

round((data.frame('ada_fdp' = ada_fdp,
                  'Bon_fdp' = bonf_fdp,
                  'cauchy_fdp' = cauchy_fdp,
                  'e_PCH_fdp' = epch_fdp,
                  'efilter_bonf_fdp' =efilter_bonf_fdp,
                  'efilter_cauchy_fdp' =efilter_cauchy_fdp)),3)

round((data.frame('ada_etp' = ada_etp,
                  'Bon_etp' = bonf_etp,
                  'cauchy_etp' = cauchy_etp,
                  'e_PCH_etp' = epch_etp,
                  'efilter_bonf_etp' = efilter_bonf_etp,
                  'efilter_cauchy_etp' =efilter_cauchy_etp)),3)

