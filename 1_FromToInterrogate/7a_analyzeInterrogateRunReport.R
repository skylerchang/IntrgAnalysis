#analyzes Interrogate run report
#expects run report in 'Data/InterrogateRunReport/' with original name
#sample IDs should have the following format:
#<SequenceId>_<OwnerLast-Patient>_<optionalGroupingVariable>
#-> only one underline; hyphens can occur more frequently -> don't use for regex search 
#18-010899-1D1P2_Piach-Mowgli_WGA+"

#Naming convention:
#pcrId, format: "15-013425-1D1P4", controls: "D16-07-1FD1P218", "NTC-can-0D1P164"
#sampleId: "k9-pc-st_D16-07.2F", k9-pc_D16-07.2F
#submissionId, format: "15-013425", controls: "D16-07", "NTC-can"

library(tidyverse)
library(here)
library(RColorBrewer)
library(reshape2)
library(gsubfn)
library(gridExtra)


#*************************** adjust settings ************************
#Is SampleNo coding? values: T/F; (e.g. MRD - yes: jake-01.1 (pbmcs), jake-01.2 (plasma); Histiocytoma - no)
sampleNoCoding<-0

#********************************************************************

targetDir<-'../../Results/InterrogateRunReport/'

files<-list.files(targetDir,pattern = '.basic--run_report.xlsx')
j<-1

for (j in 1:length(files)){
  t<-readxl::read_excel(paste0(targetDir,files[j]))
  
  #************** Run-specific code ***************
  #cleans-up inconsistencies in sample naming 
  #consider going back and fixing names in the Illumina sample sheet to make this section redundant
  #************** Run 19 - Histiocytoma ***************
  t$sample<-sub("D16-07","pc-control",t$sample)
  t$sample<-sub("NTC-can","nt-control",t$sample)
  
  
  #************** create additional sample designations ***************
  #simplify sample name -> omit '_S62_L001_R1_001'
  t$sample<-sub("_S[0-9]+_L001_R1_001","",t$sample)
  #create PcrId, format: "15-013425-1D1P4", controls: "D16-07-1FD1P218", "NTC-can-0D1P164"
  t$pcrId<-unlist(strapply(t$sample,"^(.*?)_.*$"))
  #create sampleId
  t$sampleId<-unlist(strapply(t$pcrId,"(.*)D[0-9]+P[0-9]+$"))
  #create submissionId, format: "15-013425", controls: "D16-07", "NTC-can"
  t$submissionId<-unlist(strapply(t$pcrId,"^(.*)-.*?$"))
  t$owner.patient<-unlist(strapply(t$sample,"^.*?_(.*)$"))
  #create sample fraction (MRD study: fraction1: PBMCs, fraction2: plasma; CSF study ...)
  if (sampleNoCoding){
    t$sampleFraction<-ifelse(grepl("-1D[0-9]+P[0-9]",t$pcrId),'fraction1','fraction2')
  }else{
    t$sampleFraction<-'fraction1'
  }
  t$isControl<-ifelse(grepl("control",t$submissionId),T,F)
  
  colorCount.sample<-length(unique(t$sampleId))
  getPalette<-colorRampPalette(brewer.pal(9,'Set1'))
  
  
  #*********** Run 26 specific code *************
  #this might have to be revised based on changes in sample names -> revisit
  #remove Conner's IGL samples
  t<-t[!grepl("IGL",t$sample),]
  
  #simplify sample names
  t$sample<-strapplyc(t$sample,"^(.*?_.*?)_",simplify = T)
  t$sample<-sub("-WGA_","_",t$sample)
  t$sample<-sub("-WGA","_WGA",t$sample)
  
  #***************** general stuff ******************************
  
  #clean up col names; remove white space and special characters
  colnames(t)<-gsub(" ","_",colnames(t))
  colnames(t)<-gsub("/","",colnames(t))
  colnames(t)<-gsub("[\\(\\)]","",colnames(t))
  colnames(t)<-gsub("5'","Five",colnames(t))
  colnames(t)<-gsub("3'","Three",colnames(t))
  
  
  #remove '%' from col values
  t$w_junction.percent<-as.numeric(as.character(gsub("%","",t$w_junction.percent)))
  t$Five_primed_in_R1.percent<-as.numeric(as.character(gsub("%","",t$Five_primed_in_R1.percent)))
  t$Five_primed_in_R2.percent<-as.numeric(as.character(gsub("%","",t$Five_primed_in_R2.percent)))
  t$Three_primed_in_R1.percent<-as.numeric(as.character(gsub("%","",t$Three_primed_in_R1.percent)))
  t$Three_primed_in_R2.percent<-as.numeric(as.character(gsub("%","",t$Three_primed_in_R2.percent)))
  t$w_primerdimers_in_R1.percent<-as.numeric(as.character(gsub("%","",t$w_primerdimers_in_R1.percent)))
  t$joined_normally.percent<-as.numeric(as.character(gsub("%","",t$joined_normally.percent)))
  t$joined_loosely.percent<-as.numeric(as.character(gsub("%","",t$joined_loosely.percent)))
  t$total_joined.percent<-as.numeric(as.character(gsub("%","",t$total_joined.percent)))
  t$saved_from_unjoined.percent<-as.numeric(as.character(gsub("%","",t$saved_from_unjoined.percent)))
  t$annotated.percent<-as.numeric(as.character(gsub("%","",t$annotated.percent)))
  t$w_IG_receptor.percent<-as.numeric(as.character(gsub("%","",t$w_IG_receptor.percent)))
  t$w_TR_receptor.percent<-as.numeric(as.character(gsub("%","",t$w_TR_receptor.percent)))
  t$w_A_chain.percent<-as.numeric(as.character(gsub("%","",t$w_A_chain.percent)))
  t$w_B_chain.percent<-as.numeric(as.character(gsub("%","",t$w_B_chain.percent)))
  t$w_D_chain.percent<-as.numeric(as.character(gsub("%","",t$w_D_chain.percent)))
  t$w_G_chain.percent<-as.numeric(as.character(gsub("%","",t$w_G_chain.percent)))
  t$w_H_chain.percent<-as.numeric(as.character(gsub("%","",t$w_H_chain.percent)))
  t$w_K_chain.percent<-as.numeric(as.character(gsub("%","",t$w_K_chain.percent)))
  
  #convert to numeric
  t$w_TR_receptor.count<-as.numeric(t$w_TR_receptor.count)
  t$w_A_chain.count<-as.numeric(t$w_A_chain.count)
  t$w_B_chain.count<-as.numeric(t$w_B_chain.count)
  t$w_D_chain.count<-as.numeric(t$w_D_chain.count)
  t$w_G_chain.count<-as.numeric(t$w_G_chain.count)
  t$w_H_chain.count<-as.numeric(t$w_H_chain.count)
  t$w_K_chain.count<-as.numeric(t$w_K_chain.count)
  t$diversity_normed.count<-as.numeric(t$diversity_normed.count)
  t$alpha_diversity.count<-as.numeric(t$alpha_diversity.count)
  t$effective_species.count<-as.numeric(t$effective_species.count)

  #============= plot read numbers =============
  plotlist<-list()
  x<-c('raw.count','raw.percent','w_junction.percent', 'w_junction.count')
  for (i in 1:length(x)){
    p<-ggplot(t,aes_string("sample",x[i]))+geom_boxplot()+coord_flip()
    plotlist[[i]]<-ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_readNos.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #============= plot diversity =============
  t$sampleId<-as.factor(t$sampleId)
  plotlist<-list()
  x<-c('diversity_normed.count','alpha_diversity.count','effective_species.count')
  for (i in 1:length(x)){
    #pdf(paste0(x[i],'.pdf'))
    p<-ggplot(t,aes_string("sampleId",x[i]))+geom_boxplot()+coord_flip()
    p<-p+ggtitle(x[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
    #dev.off()
  }
  pdf(paste0(targetDir,'RunReportPlots_diversity.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #============= plot priming =============
  plotlist<-list()
  x<-c('Five_primed_in_R1.percent',
       'Five_primed_in_R2.percent',
       'Three_primed_in_R1.percent',
       'Three_primed_in_R2.percent')
  titles<-c("5\' primed in R1","5\' primed in R2","3\' primed in R1","3\' primed in R2")
  for (i in 1:length(x)){
    p<-ggplot(t[t$sampleType!='control',],aes_string(x="wga",y=x[i],fill="sampleType"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ylim(0,40)+ggtitle(titles[i])+xlab('WGA')+ylab('% reads')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_priming.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #============= primer dimers=============
  plotlist<-list()
  x<-c('w_primerdimers_in_R1.percent','w_primerdimers_in_R2.percent','w_primerdimers_in_R1.count','w_primerdimers_in_R2.count')
  titles<-c("Primer dimers in R1 - percent","Primer dimers in R2 - percent","Primer dimers in R1 - count","Primer dimers in R2 - count")
  for (i in 1:length(x)){
    p<-ggplot(t[t$sampleType!='control',],aes_string(x="wga",y=x[i],fill="sampleType"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ggtitle(titles[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_primerDimers.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #============= joining ============= 
  plotlist<-list()
  x<-c("joined_normally.count","joined_normally.percent",      
        "joined_loosely.count","joined_loosely.percent",       
        "total_joined.count","total_joined.percent",      
        "saved_from_unjoined.count","saved_from_unjoined.percent")
  for (i in 1:length(x)){
    p<-ggplot(t[t$sampleType!='control',],aes_string(x="wga",y=x[i],fill="sampleType"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ggtitle(x[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_joining.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #boxplot chains
  plotlist<-list()
  x<-c("annotated.count","annotated.percent",            
       "w_IG_receptor.count","w_IG_receptor.percent",        
       "w_TR_receptor.count","w_TR_receptor.percent",         
       "w_A_chain.count","w_A_chain.percent",         
       "w_B_chain.count","w_B_chain.percent",
       "w_D_chain.count","w_D_chain.percent",
       "w_G_chain.count","w_G_chain.percent",
       "w_H_chain.count","w_H_chain.percent",
       "w_K_chain.count","w_K_chain.percent")
  for (i in 1:length(x)){
    p<-ggplot(t[t$sampleType!='control',],aes_string(x="wga",y=x[i],fill="sampleType"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ggtitle(x[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_chains.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  #=============== by locus =========================
  #create subset that includes percentage
  s<-c("sample","submission","sampleType","wga",
  "w_A_chain.percent",         
  "w_B_chain.percent",
  "w_D_chain.percent",
  "w_G_chain.percent",
  "w_H_chain.percent",
  "w_K_chain.percent")
  t.lociPercent<-t[,s]
  t.lociPercent<-melt(t.lociPercent[t.loci$sampleType!="control",],id.vars =c("sample","submission","sampleType","wga"),variable.name = "locus",value.name = 'percentage')
  
  #create subset that includes count
  s<-c("sample","submission","sampleType","wga",
       "w_A_chain.count",         
       "w_B_chain.count",
       "w_D_chain.count",
       "w_G_chain.count",
       "w_H_chain.count",
       "w_K_chain.count")
  t.lociCount<-t[,s]
  t.lociCount<-melt(t.lociCount[t.loci$sampleType!="control",],id.vars =c("sample","submission","sampleType","wga"),variable.name = "locus",value.name = 'count')
  
  
  plotlist<-list()
  #boxplot by sampleType - percent
  p<-ggplot(t.lociPercent,aes(locus,percentage,fill=sampleType))+geom_boxplot()+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))+coord_flip()+theme(legend.position="top")
  plotlist[[1]]<-ggplotGrob(p)

  #boxplot by WGA - percent
  p<-ggplot(t.lociPercent,aes(locus,percentage,fill=wga))+geom_boxplot()+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('darkolivegreen3','coral3'))+coord_flip()+theme(legend.position="top")
  plotlist[[2]]<-ggplotGrob(p)
  
  #boxplot by sampleType - count
  p<-ggplot(t.lociCount,aes(locus,count,fill=sampleType))+geom_boxplot()+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))+coord_flip()+theme(legend.position="top")
  plotlist[[3]]<-ggplotGrob(p)
  
  #boxplot by WGA - count
  p<-ggplot(t.lociCount,aes(locus,count,fill=wga))+geom_boxplot()+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('darkolivegreen3','coral3'))+coord_flip()+theme(legend.position="top")
  plotlist[[4]]<-ggplotGrob(p)
  
  pdf(paste0(targetDir,'RunReportPlots_chains_byWgaAndSampleType.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
}
i<-1
colnames(t)
