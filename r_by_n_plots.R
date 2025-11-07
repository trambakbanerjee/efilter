
library(Rfast)
library(ggplot2)
library(ggpubr)
library(latex2exp)

######################### FDR control #####################

plotter1<- function(result,rho_val,alpha_val){
  
  out = c('2/2','2/4','4/4','2/8','4/8','8/8')
  
  ada_fdp = ada_etp = bonf_fdp = bonf_etp = efilter_bonf_fdp = efilter_bonf_etp = matrix(0,6,2)
  efilter_cauchy_fdp = efilter_cauchy_etp = cauchy_fdp = cauchy_etp = epch_fdp = epch_etp = matrix(0,6,2)
  
  ada_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  ada_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  
  bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  
  efilter_bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  efilter_bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  
  efilter_cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  efilter_cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  
  cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  
  epch_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  epch_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  
  plotdata1 = data.frame('mean_fdp'=c(ada_fdp[,1],bonf_fdp[,1],
                                      cauchy_fdp[,1],efilter_bonf_fdp[,1],
                                      efilter_cauchy_fdp[,1],epch_fdp[,1]),
                         'sd_fdp'=c(ada_fdp[,2],bonf_fdp[,2],
                                    cauchy_fdp[,2],efilter_bonf_fdp[,2],
                                    efilter_cauchy_fdp[,2],epch_fdp[,2]),
                         'mean_etp'=c(ada_etp[,1],bonf_etp[,1],
                                      cauchy_etp[,1],efilter_bonf_etp[,1],
                                      efilter_cauchy_etp[,1],epch_etp[,1]),
                         'sd_etp'=c(ada_etp[,2],bonf_etp[,2],
                                    cauchy_etp[,2],efilter_bonf_etp[,2],
                                    efilter_cauchy_etp[,2],epch_etp[,2]),
                         'rbyn' = rep(out,6),
                         'method' = rep(c('AdaFilter', 'Bonferroni','Cauchy','e-Filter B',
                                          'e-Filter C', 'e-PCH'),each=6))
  
  g1<-ggplot(plotdata1, aes(x = rbyn, y = mean_fdp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = pmax(0,mean_fdp - sd_fdp), 
                      ymax = mean_fdp + sd_fdp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    geom_hline(yintercept = alpha_val,linetype = "dashed")+
    labs(title = paste0("rho = ",rho_val),
         x = "r / n",
         y = "FDR") +
    theme_minimal()
  
  g2<-ggplot(plotdata1, aes(x = rbyn, y = mean_etp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = mean_etp - sd_etp, 
                      ymax = mean_etp + sd_etp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    labs(title = paste0("rho = ",rho_val),
         x = "r / n",
         y = "Recall") +
    theme_minimal()
  
  return(list('g1'=g1,'g2'=g2))
}
plotter2<- function(result,pi0_val,alpha_val){
  
  out = c('2/2','2/4','4/4','2/8','4/8','8/8')
  
  ada_fdp = ada_etp = bonf_fdp = bonf_etp = efilter_bonf_fdp = efilter_bonf_etp = matrix(0,6,2)
  efilter_cauchy_fdp = efilter_cauchy_etp = cauchy_fdp = cauchy_etp = epch_fdp = epch_etp = matrix(0,6,2)
  
  ada_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  ada_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  
  bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  
  efilter_bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  efilter_bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  
  efilter_cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  efilter_cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  
  cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  
  epch_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  epch_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  
  plotdata1 = data.frame('mean_fdp'=c(ada_fdp[,1],bonf_fdp[,1],
                                      cauchy_fdp[,1],efilter_bonf_fdp[,1],
                                      efilter_cauchy_fdp[,1],epch_fdp[,1]),
                         'sd_fdp'=c(ada_fdp[,2],bonf_fdp[,2],
                                    cauchy_fdp[,2],efilter_bonf_fdp[,2],
                                    efilter_cauchy_fdp[,2],epch_fdp[,2]),
                         'mean_etp'=c(ada_etp[,1],bonf_etp[,1],
                                      cauchy_etp[,1],efilter_bonf_etp[,1],
                                      efilter_cauchy_etp[,1],epch_etp[,1]),
                         'sd_etp'=c(ada_etp[,2],bonf_etp[,2],
                                    cauchy_etp[,2],efilter_bonf_etp[,2],
                                    efilter_cauchy_etp[,2],epch_etp[,2]),
                         'rbyn' = rep(out,6),
                         'method' = rep(c('AdaFilter', 'Bonferroni','Cauchy','e-Filter B',
                                          'e-Filter C', 'e-PCH'),each=6))
  
  g1<-ggplot(plotdata1, aes(x = rbyn, y = mean_fdp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = pmax(0,mean_fdp - sd_fdp), 
                      ymax = mean_fdp + sd_fdp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    geom_hline(yintercept = alpha_val,linetype = "dashed")+
    labs(title = paste0("pi_00 = ",pi0_val),
         x = "r / n",
         y = "FDR") +
    theme_minimal()
  
  g2<-ggplot(plotdata1, aes(x = rbyn, y = mean_etp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = mean_etp - sd_etp, 
                      ymax = mean_etp + sd_etp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    labs(title = paste0("pi_00 = ",pi0_val),
         x = "r / n",
         y = "Recall") +
    theme_minimal()
  
  return(list('g1'=g1,'g2'=g2))
}

### Scenario 1--------------
load('scenario_1_fdr.RData')

alph = 0.2
result = list('list1'=out13_p98,'list2'=out14_p98,'list3'=out15_p98,
              'list4'=out16_p98,'list5'=out17_p98,'list6'=out18_p98)
out_rho0.2 = plotter1(result,0.2,alph)
result = list('list1'=out7_p98,'list2'=out8_p98,'list3'=out9_p98,
              'list4'=out10_p98,'list5'=out11_p98,'list6'=out12_p98)
out_rho0.4 = plotter1(result,0.4,alph)
result = list('list1'=out1_p98,'list2'=out2_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_rho0.6 = plotter1(result,0.6,alph)
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)
out_rho0.8 = plotter1(result,0.8,alph)

ggarrange(out_rho0.2$g1,out_rho0.2$g2,
          out_rho0.4$g1,out_rho0.4$g2,
          out_rho0.6$g1,out_rho0.6$g2,
          out_rho0.8$g1,out_rho0.8$g2,
          ncol=2,nrow=4,common.legend = TRUE)


### Negative Dependence --------------
load('negative_dependence_fdr.RData')

alph = 0.2
result = list('list1'=out13_p98,'list2'=out14_p98,'list3'=out15_p98,
              'list4'=out16_p98,'list5'=out17_p98,'list6'=out18_p98)
out_rho0.2 = plotter1(result,-0.2,alph)
result = list('list1'=out7_p98,'list2'=out8_p98,'list3'=out9_p98,
              'list4'=out10_p98,'list5'=out11_p98,'list6'=out12_p98)
out_rho0.4 = plotter1(result,-0.4,alph)
result = list('list1'=out1_p98,'list2'=out2_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_rho0.6 = plotter1(result,-0.6,alph)
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)
out_rho0.8 = plotter1(result,-0.8,alph)

ggarrange(out_rho0.2$g1,out_rho0.2$g2,
          out_rho0.4$g1,out_rho0.4$g2,
          out_rho0.6$g1,out_rho0.6$g2,
          out_rho0.8$g1,out_rho0.8$g2,
          ncol=2,nrow=4,common.legend = TRUE)

### Common Controls----
load('common_controls_fdr.RData')

alph = 0.2
result = list('list1'=out1_p98,'list2'=out1_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_pi098 = plotter2(result,0.98,alph)

result = list('list1'=out1_p8,'list2'=out2_p8,'list3'=out3_p8,
              'list4'=out4_p8,'list5'=out5_p8,'list6'=out6_p8)
out_pi08 = plotter2(result,0.8,alph)

ggarrange(out_pi098$g1,out_pi098$g2,
          out_pi08$g1,out_pi08$g2,
          ncol=2,nrow=2,common.legend = TRUE)

### Scale Mixture of Gaussians----
load('scale_mixture_fdr.RData')

alph = 0.2
result = list('list1'=out19_p98,'list2'=out20_p98,'list3'=out21_p98,
              'list4'=out22_p98,'list5'=out23_p98,'list6'=out24_p98)
out_rho0 = plotter1(result,0,alph)
result = list('list1'=out7_p98,'list2'=out8_p98,'list3'=out9_p98,
              'list4'=out10_p98,'list5'=out11_p98,'list6'=out12_p98)
out_rho0.4 = plotter1(result,0.4,alph)
result = list('list1'=out1_p98,'list2'=out2_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_rho0.6 = plotter1(result,0.6,alph)
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)
out_rho0.8 = plotter1(result,0.8,alph)

ggarrange(out_rho0$g1,out_rho0$g2,
          out_rho0.4$g1,out_rho0.4$g2,
          out_rho0.6$g1,out_rho0.6$g2,
          out_rho0.8$g1,out_rho0.8$g2,
          ncol=2,nrow=4,common.legend = TRUE)


##########################################################

######################### PFER control ##################
plotter1<- function(result,rho_val,alpha_val){
  
  out = c('2/2','2/4','4/4','2/8','4/8','8/8')
  
  ada_fdp = ada_etp = bonf_fdp = bonf_etp = efilter_bonf_fdp = efilter_bonf_etp = matrix(0,6,2)
  efilter_cauchy_fdp = efilter_cauchy_etp = cauchy_fdp = cauchy_etp = epch_fdp = epch_etp = matrix(0,6,2)
  
  ada_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  ada_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  
  bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  
  efilter_bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  efilter_bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  
  efilter_cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  efilter_cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  
  cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  
  epch_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  epch_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  
  plotdata1 = data.frame('mean_fdp'=c(ada_fdp[,1],bonf_fdp[,1],
                                      cauchy_fdp[,1],efilter_bonf_fdp[,1],
                                      efilter_cauchy_fdp[,1],epch_fdp[,1]),
                         'sd_fdp'=c(ada_fdp[,2],bonf_fdp[,2],
                                    cauchy_fdp[,2],efilter_bonf_fdp[,2],
                                    efilter_cauchy_fdp[,2],epch_fdp[,2]),
                         'mean_etp'=c(ada_etp[,1],bonf_etp[,1],
                                      cauchy_etp[,1],efilter_bonf_etp[,1],
                                      efilter_cauchy_etp[,1],epch_etp[,1]),
                         'sd_etp'=c(ada_etp[,2],bonf_etp[,2],
                                    cauchy_etp[,2],efilter_bonf_etp[,2],
                                    efilter_cauchy_etp[,2],epch_etp[,2]),
                         'rbyn' = rep(out,6),
                         'method' = rep(c('AdaFilter', 'Bonferroni','Cauchy','e-Filter B',
                                          'e-Filter C', 'e-PCH'),each=6))
  
  g1<-ggplot(plotdata1, aes(x = rbyn, y = mean_fdp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = pmax(0,mean_fdp - sd_fdp), 
                      ymax = mean_fdp + sd_fdp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    geom_hline(yintercept = alpha_val,linetype = "dashed")+
    labs(title = paste0("rho = ",rho_val),
         x = "r / n",
         y = "PFER") +
    theme_minimal()
  
  g2<-ggplot(plotdata1, aes(x = rbyn, y = mean_etp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = mean_etp - sd_etp, 
                      ymax = mean_etp + sd_etp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    labs(title = paste0("rho = ",rho_val),
         x = "r / n",
         y = "Recall") +
    theme_minimal()
  
  return(list('g1'=g1,'g2'=g2))
}
plotter2<- function(result,pi0_val,alpha_val){
  
  out = c('2/2','2/4','4/4','2/8','4/8','8/8')
  
  ada_fdp = ada_etp = bonf_fdp = bonf_etp = efilter_bonf_fdp = efilter_bonf_etp = matrix(0,6,2)
  efilter_cauchy_fdp = efilter_cauchy_etp = cauchy_fdp = cauchy_etp = epch_fdp = epch_etp = matrix(0,6,2)
  
  ada_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_fdp)))
  ada_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  ada_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$ada_etp)))
  
  bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_fdp)))
  bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$bpvals_etp)))
  
  efilter_bonf_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_fdp)))
  efilter_bonf_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  efilter_bonf_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_bonf_etp)))
  
  efilter_cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_fdp)))
  efilter_cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  efilter_cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$efilter_cauchy_etp)))
  
  cauchy_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_fdp)))
  cauchy_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  cauchy_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$cauchy_etp)))
  
  epch_fdp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_fdp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_fdp)))
  epch_etp[,1] = sapply(1:length(names(result)),function(k) mean(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  epch_etp[,2] = sapply(1:length(names(result)),function(k) sd(sapply(1:reps,function(i) result[[k]][[i]]$e_PCH_etp)))
  
  plotdata1 = data.frame('mean_fdp'=c(ada_fdp[,1],bonf_fdp[,1],
                                      cauchy_fdp[,1],efilter_bonf_fdp[,1],
                                      efilter_cauchy_fdp[,1],epch_fdp[,1]),
                         'sd_fdp'=c(ada_fdp[,2],bonf_fdp[,2],
                                    cauchy_fdp[,2],efilter_bonf_fdp[,2],
                                    efilter_cauchy_fdp[,2],epch_fdp[,2]),
                         'mean_etp'=c(ada_etp[,1],bonf_etp[,1],
                                      cauchy_etp[,1],efilter_bonf_etp[,1],
                                      efilter_cauchy_etp[,1],epch_etp[,1]),
                         'sd_etp'=c(ada_etp[,2],bonf_etp[,2],
                                    cauchy_etp[,2],efilter_bonf_etp[,2],
                                    efilter_cauchy_etp[,2],epch_etp[,2]),
                         'rbyn' = rep(out,6),
                         'method' = rep(c('AdaFilter', 'Bonferroni','Cauchy','e-Filter B',
                                          'e-Filter C', 'e-PCH'),each=6))
  
  g1<-ggplot(plotdata1, aes(x = rbyn, y = mean_fdp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = pmax(0,mean_fdp - sd_fdp), 
                      ymax = mean_fdp + sd_fdp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    geom_hline(yintercept = alpha_val,linetype = "dashed")+
    labs(title = paste0("pi_00 = ",pi0_val),
         x = "r / n",
         y = "PFER") +
    theme_minimal()
  
  g2<-ggplot(plotdata1, aes(x = rbyn, y = mean_etp, fill = method)) +
    geom_col(position = "dodge", color = "black") +
    geom_errorbar(aes(ymin = mean_etp - sd_etp, 
                      ymax = mean_etp + sd_etp),
                  position = position_dodge(width = 0.9), # Adjust width to match bars
                  width = 0.2) + # Width of the error bar caps
    labs(title = paste0("pi_00 = ",pi0_val),
         x = "r / n",
         y = "Recall") +
    theme_minimal()
  
  return(list('g1'=g1,'g2'=g2))
}

### Scenario 1--------------
load('scenario_1_pfer.RData')

alph = 1
result = list('list1'=out13_p98,'list2'=out14_p98,'list3'=out15_p98,
              'list4'=out16_p98,'list5'=out17_p98,'list6'=out18_p98)
out_rho0.2 = plotter1(result,0.2,alph)
result = list('list1'=out7_p98,'list2'=out8_p98,'list3'=out9_p98,
              'list4'=out10_p98,'list5'=out11_p98,'list6'=out12_p98)
out_rho0.4 = plotter1(result,0.4,alph)
result = list('list1'=out1_p98,'list2'=out2_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_rho0.6 = plotter1(result,0.6,alph)
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)
out_rho0.8 = plotter1(result,0.8,alph)

ggarrange(out_rho0.2$g1,out_rho0.2$g2,
          out_rho0.4$g1,out_rho0.4$g2,
          out_rho0.6$g1,out_rho0.6$g2,
          out_rho0.8$g1,out_rho0.8$g2,
          ncol=2,nrow=4,common.legend = TRUE)

#ggarrange(out_rho0.2$g2,out_rho0.4$g2,out_rho0.6$g2,
#          out_rho0.8$g2,ncol=2,nrow=2,common.legend = TRUE)

### Negative Dependence--------------
load('negative_dependence_pfer.RData')

alph = 1
result = list('list1'=out13_p98,'list2'=out14_p98,'list3'=out15_p98,
              'list4'=out16_p98,'list5'=out17_p98,'list6'=out18_p98)
out_rho0.2 = plotter1(result,-0.2,alph)
result = list('list1'=out7_p98,'list2'=out8_p98,'list3'=out9_p98,
              'list4'=out10_p98,'list5'=out11_p98,'list6'=out12_p98)
out_rho0.4 = plotter1(result,-0.4,alph)
result = list('list1'=out1_p98,'list2'=out2_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_rho0.6 = plotter1(result,-0.6,alph)
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)
out_rho0.8 = plotter1(result,-0.8,alph)

ggarrange(out_rho0.2$g1,out_rho0.2$g2,
          out_rho0.4$g1,out_rho0.4$g2,
          out_rho0.6$g1,out_rho0.6$g2,
          out_rho0.8$g1,out_rho0.8$g2,
          ncol=2,nrow=4,common.legend = TRUE)

### Common Controls----
load('common_controls_pfer.RData')

alph = 1
result = list('list1'=out1_p98,'list2'=out1_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_pi098 = plotter2(result,0.98,alph)

result = list('list1'=out1_p8,'list2'=out2_p8,'list3'=out3_p8,
              'list4'=out4_p8,'list5'=out5_p8,'list6'=out6_p8)
out_pi08 = plotter2(result,0.8,alph)

ggarrange(out_pi098$g1,out_pi098$g2,
          out_pi08$g1,out_pi08$g2,
          ncol=2,nrow=2,common.legend = TRUE)

### Scale Mixture of Gaussians----
load('scale_mixture_pfer.RData')

alph = 1
result = list('list1'=out19_p98,'list2'=out20_p98,'list3'=out21_p98,
              'list4'=out22_p98,'list5'=out23_p98,'list6'=out24_p98)
out_rho0 = plotter1(result,0,alph)
result = list('list1'=out7_p98,'list2'=out8_p98,'list3'=out9_p98,
              'list4'=out10_p98,'list5'=out11_p98,'list6'=out12_p98)
out_rho0.4 = plotter1(result,0.4,alph)
result = list('list1'=out1_p98,'list2'=out2_p98,'list3'=out3_p98,
              'list4'=out4_p98,'list5'=out5_p98,'list6'=out6_p98)
out_rho0.6 = plotter1(result,0.6,alph)
result = list('list1'=outa_p98,'list2'=outb_p98,'list3'=outc_p98,
              'list4'=outd_p98,'list5'=oute_p98,'list6'=outf_p98)
out_rho0.8 = plotter1(result,0.8,alph)

ggarrange(out_rho0$g1,out_rho0$g2,
          out_rho0.4$g1,out_rho0.4$g2,
          out_rho0.6$g1,out_rho0.6$g2,
          out_rho0.8$g1,out_rho0.8$g2,
          ncol=2,nrow=4,common.legend = TRUE)
