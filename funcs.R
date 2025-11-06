

library(Rfast)
library(mvtnorm)

e_filter<-function(eS,eF,m,alph,TypeI){
  
  idx = order(eS,decreasing = TRUE)
  eS_ord = eS[idx]
  #m_ord = sapply(1:m,function(i) sum(1*(eF>=eS_ord[i])))
  
  ef_sorted <- sort(eF)
  m_ord = m - findInterval(eS_ord, ef_sorted, left.open = TRUE)
  if(TypeI=="FDR"){
    
    temp = eS_ord*(1:m)/m_ord
    #e_ord = matrix(0,m,1)
    #e_ord[1] = max(temp)
    #e_ord[2:m] = sapply(2:m,function(i) max(temp[i:m]))
    e_ord = rev(cummax(rev(temp)))
    
  }
  if(TypeI=="PFER"){
    
    e_ord = eS_ord/m_ord
  }
  de_ord = (e_ord>(1/alph))
  de = evals = matrix(0,m,1)
  de[idx] = 1*de_ord
  evals[idx] = e_ord
  
  return(list("decisions"=de,'evals'=evals))
  
}

e_PCH<-function(evals,m,n,r,alph,TypeI){
  
  pch_evals = matrix(0,m,1)
  for(i in 1:m){
    
    cc = evals[i,]
    cc = cc[order(cc)]
    pch_evals[i] = (1/(n-r+1))*sum(cc[1:(n-r+1)])
  }
  if(TypeI=="FDR"){
    
    ebh = ebh.func(pch_evals,alph)
    de = ebh$de
  }
  if(TypeI=="PFER"){
    
    de = 1*(pch_evals>=(m/alph))
  }
  
  return(list("decisions"=de,'evals'=evals))
}

cauchy_PCH<-function(pvals,m,n,r,alph,TypeI){
  
  cauchypvals = cauchy_eF = matrix(0,m,1)
  for(i in 1:m){
    
    cc = pvals[i,]
    cc = cc[order(cc)]
    test_stat = (1/(n-r+1))*sum(tan((0.5-cc[r:n])*pi))
    cauchypvals[i] = pcauchy(test_stat,lower.tail = FALSE)
    test_stat_1 = tan((0.5-cc[r-1])*pi)#(1/(agents-r+1))*sum(tan((0.5-cc[(r-1):agents])*pi))
    cauchy_eF[i] = pcauchy(test_stat_1,lower.tail = FALSE)
  }
  if(TypeI=="FDR"){
    
    de = bh.func(cauchypvals,alph)$de
  }
  if(TypeI=="PFER"){
    
    de = 1*(cauchypvals<=(alph/m))
  }
  
  return(list("decisions"=de,'pvals'=cauchypvals,'cauchy_eF'=cauchy_eF))
}

ebh.func<-function(ev, q)
{ 
  # the input: 
  # ev: the e-values
  # q: the control FDR level
  # the output :
  # thr: the e-value threshold
  # de: the decision rule
  
  m=length(ev)
  st.ev<-sort(ev,decreasing = TRUE)   
  evi<-st.ev*1:m
  hps<-rep(0,m)
  k<-max(which(evi>=(m/q)))
  ek<-st.ev[k]
  hps[which(ev>=ek)]<-1
  
  return (list(thr=ek, de=hps))
}

bh.func<-function(pv, q)
{ 
  # the input 
  # pv: the p-values
  # q: the FDR level
  # the output 
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # re: the index of rejected hypotheses
  # ac: the index of accepted hypotheses
  # de: the decision rule
  
  m=length(pv)
  st.pv<-sort(pv)   
  #print(length(st.pv))
  #print(m)
  pvi<-st.pv/1:m
  hps<-rep(0, m)
  if (max(pvi<=(q/m))==0)
  {
    k<-0
    pk<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=pk, re=reject, ac=accept, de=hps)
  return (y)
}

GenPMat_corr <- function(M = 1000,Sigma_chol.1,Sigma_chol.2,
                         n = 2,
                         r = 2,
                         all.zero.frac = 0.8,
                         alternative.frac = 0.1,
                         mu = c(3,-3,4,-4,6,-6),
                         one.sided = F) 
{
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  
  noise.matrix = matrix(rnorm(M*n),M,n)
  
  ## generate dependency across screening tests and agents
  U = Sigma_chol.1#cholesky(rho1*matrix(1,M,M)+(1-rho1)*diag(M))
  V = Sigma_chol.2#chol(rho2*matrix(1,n,n)+(1-rho2)*diag(n))
  
  zmat = mean.matrix+t(U)%*%noise.matrix%*%V
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}

GenPMat_corr_scalemix <- function(M = 1000,Sigma_chol.1,Sigma_chol.2,
                         n = 2,
                         r = 2,
                         all.zero.frac = 0.8,
                         alternative.frac = 0.1,
                         mu = c(3,-3,4,-4,6,-6),
                         one.sided = F) 
{
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  cc = c(rnorm(ceiling(M*n/2),0,sqrt(1)),
         rnorm(ceiling(M*n/4),0,sqrt(2)),
         rnorm(M*n-3*ceiling(M*n/4),0,sqrt(4)))
  
  noise.matrix = matrix(cc,M,n)
  
  ## generate dependency across screening tests and agents
  U = Sigma_chol.1#cholesky(rho1*matrix(1,M,M)+(1-rho1)*diag(M))
  V = Sigma_chol.2#chol(rho2*matrix(1,n,n)+(1-rho2)*diag(n))
  
  zmat = mean.matrix+t(U)%*%noise.matrix%*%V
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}

GenPMat_corr_t <- function(M = 1000,Sigma,
                           n = 2,
                           r = 2,
                           t.df = 5,
                           all.zero.frac = 0.8,
                           alternative.frac = 0.1,
                           mu = c(3,-3,4,-4,6,-6),
                           one.sided = F) 
{
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  
  noise.matrix = mvtnorm::rmvt(M,Sigma,df=t.df,type="Kshirsagar")#t_columns_with_corr(M,n,t.df,Sigma_chol)#matrix(rnorm(M*n),M,n)
  
  zmat = mean.matrix+noise.matrix
  
  # signs <- 2 * (zmat >= 0) - 1
  # raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  # 
  # if (one.sided) {
  #   pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  # } else {
  #   pvalue.matrix <- 2 * raw.pvalues
  # }
 pvalue.matrix = 2*(1-pt(abs(zmat),t.df,ncp=0))
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}


GenPMat_common_control <- function(M = 1000,
                         n = 2,
                         r = 2,
                         all.zero.frac = 0.8,
                         alternative.frac = 0.01,
                         mu = c(3,-3,4,-4,6,-6),
                         one.sided = F) 
{
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  
  control.y = matrix(rnorm(M),M,1)
  treat.x = matrix(0,M,n)
  zmat = matrix(0,M,n)
  for(ii in 1:n){
    set.seed(ii)
    treat.x[,ii] = sapply(1:M,function(i) rnorm(1,mean.matrix[i,ii],1))
    zmat[,ii] = (treat.x[,ii]-control.y)/sqrt(2)
  }
  ##noise.matrix = matrix(rnorm(M*n),M,n)
  
  ## generate dependency across screening tests and agents
  ##U = Sigma_chol.1#cholesky(rho1*matrix(1,M,M)+(1-rho1)*diag(M))
  ##V = Sigma_chol.2#chol(rho2*matrix(1,n,n)+(1-rho2)*diag(n))
  
  ##zmat = #mean.matrix+t(U)%*%noise.matrix%*%V
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}

GenPMat_commoncases_n4 <- function(M = 1000,
                                   r = 2,
                                   s = 10,
                                   all.zero.frac = 0.8,
                                   alternative.frac = 0.01,
                                   mu = c(3,-3,4,-4,6,-6),
                                   one.sided = F) 
{
  n = 4
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  noise.mat = matrix(0,M,n)
  set.seed(1)
  temp1 = matrix(rnorm(s*M,0,sqrt(s)),s,M)
  noise.mat[,1] = colmeans(temp1)
  set.seed(2)
  temp2 = matrix(rnorm(s*M,0,sqrt(s)),s,M)
  noise.mat[,2] = colmeans(rbind(temp1,temp2))
  set.seed(3)
  temp3 = matrix(rnorm(s*M,0,sqrt(s)),s,M)
  noise.mat[,3] = colmeans(rbind(temp3,
                                 temp2))
  set.seed(4)
  noise.mat[,4] = colmeans(rbind(temp1,
                                 temp2,
                                 temp3))
  
  zmat = mean.matrix+noise.mat
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}

GenPMat_commoncases_n5 <- function(M = 1000,
                                   r = 2,
                                   s = 5000,
                                   q = 400,
                                   all.zero.frac = 0.8,
                                   alternative.frac = 0.01,
                                   mu = c(3,-3,4,-4,6,-6),
                                   one.sided = F) 
{
  n = 5
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  noise.mat = matrix(0,M,n)
  temp1 = matrix(rnorm(s*M,0,sqrt(s)),s,M)
  noise.mat[,1] = colmeans(temp1) # Study 1
  temp2 = sqrt(1000)*temp1[1:1000,]/sqrt(s)# Study 2
  noise.mat[,2] = colmeans(temp2)
  temp3 = sqrt(1000)*temp1[1001:2000,]/sqrt(s)# Study 3 
  noise.mat[,3] = colmeans(temp3) 
  temp4 = sqrt(1000)*temp1[2001:3000,]/sqrt(s)# Study 4
  noise.mat[,4] = colmeans(temp4)
  temp5 = rbind(temp2[1:q,],
                matrix(rnorm((1000-q)*M,0,sqrt(1000)),(1000-q),M)) # Study 5
  noise.mat[,5] = colmeans(temp5)
  
  zmat = mean.matrix+noise.mat
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}

GenPMat_commoncases_n2 <- function(M = 1000,
                                   r = 2,
                                   s = 1000,
                                   all.zero.frac = 0.8,
                                   alternative.frac = 0.01,
                                   mu = c(3,-3,4,-4,6,-6),
                                   one.sided = F) 
{
  n = 2
  combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
                                                       n), nrow = 2))))
  category <- rowSums(combs)
  n.cases <- sum(category >= r)
  weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
                 nrow(combs))
  weights[category == 0] <- all.zero.frac
  weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)
  
  mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]
  
  
  
  truth.pc <- rowSums(mean.matrix) >= r
  
  
  ## This is the signal matrix
  mean.matrix <- apply(mean.matrix, 1, function(v) {
    
    v[v == 1] <- sample(mu, sum(v == 1), replace = T)
    return(v)
  })
  mean.matrix <- t(mean.matrix)
  
  noise.mat = matrix(0,M,n)
  set.seed(1)
  temp1 = matrix(rnorm(s*M,0,sqrt(s)),s,M)
  noise.mat[,1] = colmeans(temp1)
  set.seed(2)
  temp2 = matrix(rnorm(50*M,0,sqrt(s)),s,M)
  noise.mat[,2] = colmeans(rbind(temp1[1:950,],temp2))
  # set.seed(3)
  # temp3 = matrix(rnorm(s*M,0,sqrt(s)),s,M)
  # noise.mat[,3] = colmeans(rbind(temp3,
  #                                temp2))
  # set.seed(4)
  # noise.mat[,4] = colmeans(rbind(temp1,
  #                                temp2,
  #                                temp3))
  
  zmat = mean.matrix+noise.mat
  
  signs <- 2 * (zmat >= 0) - 1
  raw.pvalues <- pnorm(abs(zmat), lower.tail = F)
  
  if (one.sided) {
    pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
  } else {
    pvalue.matrix <- 2 * raw.pvalues 
  }
  
  return(list(truth.pc = truth.pc, 
              pvalue.mat = pvalue.matrix,
              zmat = zmat))
}

AR1<-function(p,rho=0.5)
{
  Sigma<-matrix(0,p,p);
  for (i in 1:p)
  {
    for (j in 1:p)
    {
      Sigma[i,j]<-rho^(abs(i-j));
    }
  }
  return(Sigma);
}

AR1_chol<-function(p,rho=0.5)
{
  R = matrix(0,p,p);                # allocate p x p matrix
  R[1,] = rho^(0:(p-1));        # formula for 1st row
  cc = sqrt(1 - rho^2); # scaling factor: c^2 + rho^2 = 1
  R2 = cc*R[1,];        # formula for 2nd row
  for(j in 2 : p){    # shift elements in 2nd row for remaining rows
    R[j, j:p] = R2[1:(p-j+1)]
  }
  return(R)
}

