#analyzes Interrogate run report
#expects run report in 'Data/InterrogateRunReport/' with original name

library(tidyverse)
library(here)
library(RColorBrewer)
library(gridExtra)
library(reshape2)

targetDir<-'../../Data/InterrogateRunReport/'

files<-list.files(targetDir,pattern = '.basic--run_report.xlsx')


for (j in 1:length(files)){
  t<-readxl::read_excel(paste0(targetDir,files[j]))
  
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
  
  
  
  
  #define additional columns
  t$sample<-gsub("_S[0-9]*_L001_R1_001","",t$sample)
  t$owner.patient<-gsub("^.*_","",t$sample)
  t$submission<-gsub("_.*$","",t$sample)
  
  t$sampleType<-ifelse(grepl("-1D[0-9]+P[0-9]+C[0-9]+L[0-9]+P[0-9]+",t$submission),'pellet','supernatant')
  t$sampleType[grepl("k9",t$submission)]<-'control'
  t$wga<-ifelse(grepl("WGA",t$submission),TRUE,FALSE)
  
  
  ggplot(t,aes(sample,raw.count))+geom_col()+coord_flip()
  ggplot(t,aes(sample,w_junction.percent,fill=sampleType))+geom_col()+coord_flip()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  
  
  #boxplots diversity
  plotlist<-list()
  x<-c('w_junction.percent', 'diversity_normed.count','alpha_diversity.count','effective_species.count')
  for (i in 1:length(x)){
    #pdf(paste0(x[i],'.pdf'))
    p<-ggplot(t[t$sampleType!='control',],aes_string(x="wga",y=x[i],fill="sampleType"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ggtitle(x[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
    #dev.off()
  }
  pdf(paste0(targetDir,'RunReportPlots_diversity.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #boxplot by priming
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
  
  #boxplot primer dimers
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
  
  #boxplot joining
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
