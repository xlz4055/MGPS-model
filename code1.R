

setwd("~/Desktop/model")
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)

iddf=read.table("imc731id.txt",header = T,sep = "\t")

imc0=as.vector(iddf$id)

result=data.frame()

imc=imc0[1:731]  
for (i in imc ) {
  expo_rt<-extract_instruments(outcome=i,p1 = 5e-6,
                               clump = T,
                               p2 = 5e-6,
                               r2 = 0.001,
                               kb = 10000) 
  
  outc_rt <- extract_outcome_data(
    snps = expo_rt$SNP,
    outcomes = '#####'
  ) 
  
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  if(result_or$pval[3]<0.05){
    result=rbind(result,cbind(id=i,pvalue=result_or$pval[3]))
    
  }
}



#######figure----
library(TwoSampleMR)
library(ggplot2)

setwd("~/Desktop/model")

iddf=read.table("result.txt",header = T,sep = "\t")   #有19个细胞

imc0=as.vector(iddf$id)
imc=imc0[1:731]

for (i in imc ) {
  expo_rt<-extract_instruments(outcome=i,p1 = 5e-6,
                               clump = T,
                               p2 = 5e-6,
                               r2 = 0.001,
                               kb = 10000) 
  
  outc_rt <- extract_outcome_data(
    snps = expo_rt$SNP,
    outcomes = '#######'
  ) 
  
  #####################################
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  dir.create(i) 
  filename=i
  
  write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
  write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
  pleiotropy=mr_pleiotropy_test(harm_rt)
  write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
  heterogeneity=mr_heterogeneity(harm_rt)
  write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
  presso=run_mr_presso(harm_rt,NbDistribution = 1000)
  capture.output(presso,file = paste0(filename,"/presso.txt"))
  
  ######  figure###############################
  p1 <- mr_scatter_plot(mr_result, harm_rt)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
  
  singlesnp_res<- mr_singlesnp(harm_rt)
  singlesnpOR=generate_odds_ratios(singlesnp_res)
  
  write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
  sen_res<- mr_leaveoneout(harm_rt)
  p3 <- mr_leaveoneout_plot(sen_res)
  
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
  res_single <- mr_singlesnp(harm_rt)
  p4 <- mr_funnel_plot(singlesnp_res)
  
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
}



















