
library(Rfast)
library(foreach)
library(doParallel)

source('funcs.R')

adafilter_x<-function (p.matrix, r, type.I.err = c("FDR", "FWER", "PFER"), 
                      alpha = 0.05, fast = F) 
{
  type.I.err <- match.arg(type.I.err, c("FDR", "FWER", "PFER"))
  n <- ncol(p.matrix)
  n.noNA <- rowSums(!is.na(p.matrix))
  if (fast) {
    sh.idx <- which((n.noNA >= r) & rowSums(p.matrix <= alpha/(n.noNA - 
                                                                 r + 1), na.rm = T) >= (r - 1))
  }
  else sh.idx <- which(n.noNA >= r)
  if (length(sh.idx) == 0) {
    if (fast) {
      warning("All p-values are too large! No hypotheses can be rejected. Return NA")
    }
    else warning(paste("No hypotheses have at least", r, 
                       "non-missing values. return NA"))
    return(NA)
  }
  else {
    p.mat.sh <- p.matrix[sh.idx, ]
    M <- sum(n.noNA >= r)
    N <- length(sh.idx)
    n.noNA <- n.noNA[sh.idx]
    p.mat.sh <- t(apply(p.mat.sh, 1, sort, method = "quick", 
                        na.last = T))
    selection.p <- p.mat.sh[, r] * (n.noNA - r + 1)
    filter.p <- p.mat.sh[, r - 1] * (n.noNA - r + 1)
    filter.p[(n.noNA < r) & (n.noNA >= r - 1)] <- NA
    sorted.s <- sort(selection.p, method = "quick", index.return = T)
    names(sorted.s$x) <- paste0("S", 1:N)
    names(filter.p) <- paste0("F", 1:N)
    temp <- sort(c(filter.p, sorted.s$x), method = "quick")
    temp[1:(2 * N)] <- 1:(2 * N)
    adj.number <- temp[names(sorted.s$x)] - 1:N
    if (type.I.err != "FDR") {
      adjusted.p <- adj.number * sorted.s$x
      adjusted.p <- cummin(adjusted.p[N:1])[N:1]
      if (type.I.err == "FWER") 
        adjusted.p <- pmin(1, adjusted.p)
    }
    else {
      adjusted.p <- adj.number/(1:N) * sorted.s$x
      adjusted.p <- pmin(1, cummin(adjusted.p[N:1])[N:1])
    }
    adj.number[sorted.s$ix] <- adj.number
    adjusted.p[sorted.s$ix] <- adjusted.p
    decision <- adjusted.p <= alpha
    results <- data.frame(matrix(NA, nrow(p.matrix), 5))
    colnames(results) <- c("decision", "adjusted.p", "selection.p", 
                           "filter.p", "adj.number")
    results[sh.idx, ] <- cbind(decision, adjusted.p, selection.p, 
                               filter.p, adj.number)
    #rownames(results) <- rownames(p.matrix)
    return(results)
  }
}

merged_data <- readRDS("merged_data.rds")
pval_mat = as.matrix(merged_data[,c('pvalue','pval1','pval2',
                                    'pval3','pval4')])

#2. Begin analysis

analysis = function(pval_mat,alpha,r){
  
  m = dim(pval_mat)[1]
  agents = dim(pval_mat)[2]
  
  # AdaFilter
  adaFilter_result <- adafilter_x(pval_mat, r = r,alpha = alpha,type.I.err = "FDR")
  
  ## Bonferroni pvals
  bpvals = matrix(0,m,1)
  F_j = matrix(0,m,1)
  for(i in 1:m){
    
    cc = c(pval_mat[i,])
    cc = cc[order(cc)]
    bpvals[i] = min(1,(agents-r+1)*cc[r])
    F_j[i] = min(1,(agents-r+1)*cc[r-1])
  }
  bh_bpvals = bh.func(bpvals,alpha)
  
  ## E-AdaFilter
  kappa_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7)
  e_ada_rejs = matrix(0,length(kappa_vals))
  for(k in 1:length(kappa_vals)){
    
    kappa = kappa_vals[k]
    e1 = kappa*(F_j^{kappa-1}) 
    e2 = kappa*(bpvals^{kappa-1}) 
    
    ff<-function(gam,e1,e2,alpha){
      
      return(gam*sum(1*(e1>(1/gam)))/max(1,sum(1*(e2>(1/gam))))-alpha)
    }
    
    sol = try(uniroot(ff,c(0,alpha),e1=e1,e2=e2,alpha=alpha))
    gam_BH = if(inherits(sol,"try-error")) NA else sol$root
    #gam_BH = uniroot(ff,c(0,alph),e1=e1,e2=e2,alpha=alph)$root
    
    
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
    
    return(gam*sum(1*(e1>(1/gam)))/max(1,sum(1*(e2>(1/gam))))-alpha)
  }
  gam_BH = uniroot(ff,c(0,alpha),e1=e1,e2=e2,alpha=alpha)$root
  
  de = 1*(e2>(1/gam_BH))
  
  ## Naive
  pch_evals_rejs = matrix(0,length(kappa_vals),1)
  for(k in 1:length(kappa_vals)){
    
    kappa = kappa_vals[k]
    evals = kappa*(pval_mat^{kappa-1}) 
    pch_evals = matrix(0,m,1)
    
    for(i in 1:m){
      
      cc = evals[i,]
      cc = cc[order(cc)]
      pch_evals[i] = (1/(agents-r+1))*sum(cc[1:(agents-r+1)])#prod(cc[1:(agents-r+1)])#
    }
    pch_evals_rejs[k] = sum(ebh.func(pch_evals,alpha)$de)
  }
  kappa = kappa_vals[which.max(pch_evals_rejs)]
  evals = kappa*(pval_mat^{kappa-1}) 
  pch_evals = matrix(0,m,1)
  for(i in 1:m){
    
    cc = evals[i,]
    cc = cc[order(cc)]
    pch_evals[i] = (1/(agents-r+1))*sum(cc[1:(agents-r+1)])#prod(cc[1:(agents-r+1)])#
  }
  ebh = ebh.func(pch_evals,alpha)
  
  
  ###############################################################
  
  return(list('rejects_eadafilter' = de,
              'rejects_enaive' = ebh$de,
              'rejects_bpvals' = bh_bpvals$de,
         'rejects_adaFilter' = adaFilter_result$decision))
}

#-------------------------------------------------------------

output_r2 = analysis(pval_mat,alpha = 0.01, r = 2)
output_r3 = analysis(pval_mat,alpha = 0.01, r = 3)
output_r4 = analysis(pval_mat,alpha = 0.01, r = 4)
output_r5 = analysis(pval_mat,alpha = 0.01, r = 5)

save.image('rep_gwas.RData')


##############################################################################

load('rep_gwas.RData')

result = output_r2
fuma_adafilter = data.frame('rsID'=merged_data$rsID,'pval'=merged_data$METAL_Pvalue,
                            'decisions'=result$rejects_adaFilter)
 
fuma_adafilter = fuma_adafilter[fuma_adafilter$decisions==1,'rsID']

write.table(fuma_adafilter, file = "adafilter_rsid_r2.txt",row.names = FALSE,
            quote=FALSE, sep = "\t")

fuma_efilter = data.frame('rsID'=merged_data$rsID,'pval'=merged_data$METAL_Pvalue,
                            'decisions'=result$rejects_eadafilter)

fuma_efilter = fuma_efilter[fuma_efilter$decisions==1,'rsID']

write.table(fuma_efilter, file = "efilter_rsid_r2.txt",row.names = FALSE,
            quote=FALSE,sep = "\t")

fuma_bonf = data.frame('rsID'=merged_data$rsID,'pval'=merged_data$METAL_Pvalue,
                          'decisions'=result$rejects_bpvals)

fuma_bonf = fuma_bonf[fuma_bonf$decisions==1,'rsID']

write.table(fuma_bonf, file = "bonf_rsid_r2.txt",row.names = FALSE,
            quote=FALSE,sep = "\t")

########################### Reading mapped genes ##################

#### r = 2
library(readr)
adafilter_r2_mappedgenes <- read_delim("adafilter_r2_mappedgenes_grch37.txt", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

adafilter_r2_mappedgenes = adafilter_r2_mappedgenes$SYMBOL
idx = which(adafilter_r2_mappedgenes=="-")
adafilter_r2_mappedgenes = adafilter_r2_mappedgenes[-idx]
adafilter_r2_mappedgenes = unique(adafilter_r2_mappedgenes)

write.table(adafilter_r2_mappedgenes, file = "adafilter_r2_mappedgenes_clean_grch37.txt",row.names = FALSE,
            quote=FALSE,sep = "\t")

efilter_r2_mappedgenes <- read_delim("efilter_r2_mappedgenes_grch37.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

efilter_r2_mappedgenes = efilter_r2_mappedgenes$SYMBOL
idx = which(efilter_r2_mappedgenes=="-")
efilter_r2_mappedgenes = efilter_r2_mappedgenes[-idx]
efilter_r2_mappedgenes = unique(efilter_r2_mappedgenes)
write.table(efilter_r2_mappedgenes, file = "efilter_r2_mappedgenes_clean_grch37.txt",row.names = FALSE,
            quote=FALSE,sep = "\t")

bonf_r2_mappedgenes <- read_delim("bonf_r2_mappedgenes_grch37.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

bonf_r2_mappedgenes = bonf_r2_mappedgenes$SYMBOL
idx = which(bonf_r2_mappedgenes=="-")
bonf_r2_mappedgenes = bonf_r2_mappedgenes[-idx]
bonf_r2_mappedgenes = unique(bonf_r2_mappedgenes)
write.table(bonf_r2_mappedgenes, file = "bonf_r2_mappedgenes_clean_grch37.txt",row.names = FALSE,
            quote=FALSE,sep = "\t")

