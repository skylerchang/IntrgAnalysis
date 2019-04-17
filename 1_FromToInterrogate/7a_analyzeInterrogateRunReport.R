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
sampleNoCoding<-1

#********************************************************************
setwd(here())
getwd()

targetDir<-'../Results/InterrogateRunReport/'

files<-list.files(targetDir,pattern = "^Run.*xlsx")
#pick which file to use
# 1 = csf run 2 = seq run 24 
# 2 = csf run 3 = seq run 26 
# 3 = excel file with both csf run 2 and 3 combined 
j<-1

for (j in 1:length(files)){
  t<-readxl::read_excel(paste0(targetDir,files[j]))

  #************** Run-specific code ***************
  #cleans-up inconsistencies in sample naming 
  #consider going back and fixing names in the Illumina sample sheet to make this section redundant
  #************** Run 19 - Histiocytoma ***************
  t$sample<-sub("D16-07","pc-control",t$sample)
  t$sample<-sub("NTC-can","nt-control",t$sample)
  #************** Run 26 - CSF ***************
  #remove controls - quick & dirty fix until run has been re-processed with new sample sheet
  t<-t[!grepl("k9",t$sample),]
  t<-t[!grepl("IGL",t$sample),]
  t$wga<-ifelse(grepl("WGA",t$sample),T,F)
  
  #************** create additional sample designations ***************
  #simplify sample name -> omit '_S62_L001_R1_001'
  t$sample<-sub("_S[0-9]+_L001_R1_001","",t$sample)
  #create PcrId, format: "15-013425-1D1P4", controls: "D16-07-1FD1P218", "NTC-can-0D1P164"
  t$pcrId<-unlist(strapply(t$sample,"^(.*?)_.*$"))
  #create sampleId
  t$sampleId<-unlist(strapply(t$pcrId,"(.*)D[0-9].*P[0-9].*$"))
  #create submissionId, format: "15-013425", controls: "D16-07", "NTC-can"
  t$submissionId<-unlist(strapply(t$sampleId,"^(.*)-.*?$"))
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
  t$raw.percent<-as.numeric(as.character(gsub("%","",t$raw.percent)))
  t$w_junction.percent<-as.numeric(as.character(gsub("%","",t$w_junction.percent)))
  t$usable.percent<-as.numeric(as.character(gsub("%","",t$usable.percent)))
  t$w_most_pop_clns.percent<-as.numeric(as.character(gsub("%","",t$w_most_pop_clns.percent)))
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
  t$joined_loosely.count<-as.numeric(t$joined_loosely.count)
  t$saved_from_unjoined.count<-as.numeric(t$saved_from_unjoined.count)
  t$raw.percent<-as.numeric(t$raw.percent)
  t$usable.percent<-as.numeric(t$usable.percent)
  t$w_most_pop_clns.percent<-as.numeric(t$w_most_pop_clns.percent)
  colnames(t)
}
  #============= plot read numbers =============
  plotlist<-list()
  x<-c('raw.count','raw.percent','w_junction.percent', 'w_junction.count')
  for (i in 1:length(x)){
    p<-ggplot(t,aes_string("sampleId",x[i]))+geom_boxplot(aes(fill=wga))+coord_flip()
    plotlist[[i]]<-ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_readNos.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #============= plot diversity - by sampleId =============
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
  pdf(paste0(targetDir,'RunReportPlots_diversity-bySampleId.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #============= plot diversity - by WGA =============
  plotlist<-list()
  x<-c('diversity_normed.count','alpha_diversity.count','effective_species.count')
  titles<-c('diversity_normed.count','alpha_diversity.count','effective_species.count')
  for (i in 1:length(x)){
    p<-ggplot(t[t$sampleFraction!='control',],aes_string(x="wga",y=x[i],fill="sampleFraction"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ylim(0,40)+ggtitle(titles[i])+xlab('WGA')+ylab('% reads')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_diversity-byWGA.pdf'))
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
    p<-ggplot(t[t$sampleFraction!='control',],aes_string(x="wga",y=x[i],fill="sampleFraction"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
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
    p<-ggplot(t[t$sampleFraction!='control',],aes_string(x="wga",y=x[i],fill="sampleFraction"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
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
    p<-ggplot(t[t$sampleFraction!='control',],aes_string(x="wga",y=x[i],fill="sampleFraction"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ggtitle(x[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_joining.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  
  #=============== by locus - count =========================
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
    p<-ggplot(t[t$sampleFraction!='control',],aes_string(x="wga",y=x[i],fill="sampleFraction"))+geom_boxplot(outlier.size = 0)+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.2))+scale_fill_manual(values=c('gray80','dodgerblue3'))
    p<-p+ggtitle(x[i])+xlab('WGA')+theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- ggplotGrob(p)
  }
  pdf(paste0(targetDir,'RunReportPlots_chains.pdf'))
  marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
  #=============== by locus - percent =========================
  #create subset that includes percentage
  s<-c("sample","submission","sampleType","wga",
  "w_A_chain.percent",         
  "w_B_chain.percent",
  "w_D_chain.percent",
  "w_G_chain.percent",
  "w_H_chain.percent",
  "w_K_chain.percent")
  t.lociPercent<-t[,s]
  t.lociPercent<-melt(t.lociPercent[t.loci$sampleFraction!="control",],id.vars =c("sample","submission","sampleType","wga"),variable.name = "locus",value.name = 'percentage')
  
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
#}
i<-1
colnames(t)



#====================== WGA+ vs WGA- ===== (word from Tamara: i am a amatuer at R so i made this very simple)
#to compare the correlations of variable and significance (t-test) data subsets will be made for WGA 
# always double check on which file "seq run" j is 

#========= subsets for ca and cf fractions for wga and non-wga 
#to trun fraction in to a numeric value -> fraction1 to 1 
#create sample fraction (MRD study: fraction1: PBMCs, fraction2: plasma; CSF study ...)
#create sample fraction (MRD study: fraction1: PBMCs, fraction2: plasma; CSF study ...)
if (sampleNoCoding){
  t$sampleFraction<-ifelse(grepl("-1D[0-9]+P[0-9]",t$pcrId),'1','2')
}else{
  t$sampleFraction<-'1'
}

# creating subset -> variable wga was already created above 
T.wga <-t[t$wga,]
F.wga <-t[!t$wga,]

#csf fraction for wga
T.wga.ca<- subset(T.wga,sampleFraction==1)
T.wga.cf <-subset(T.wga,sampleFraction==2)

#csf fraction for non-wga
F.wga.ca<- subset(F.wga,sampleFraction==1)
F.wga.cf <-subset(F.wga,sampleFraction==2)

#=============correlation overall / correlation matrix
# removing columns not of intrest to this intial data analysis 
#allows for the creation of a correlation matrix 
#to get this to work you have to remove all % signs from the varibales 
#tcor<-t[,-c(1:6,8,11,12,14,20:83,85:90)]
#plot(tcor)
#cor(tcor) 
    #something is numeric = not working 
    #above code only works for run 26 (run 24 has extra columns)

##subset for wga corelation -> just looking at the above variables 
#T.wga <-tcor[tcor$wga,]
#F.wga <-tcor[!tcor$wga,]

#correlations, means, plots and T-test(to test the means)
    #basic code to be used-> will repeat this for the multiple varibales to be looked at

#=========================================================================
#====================== raw reads - count ================================
#=========================================================================

#raw counts vs effective species
#overall
cor.test(t$raw.count,t$effective_species.count)
plot(t$raw.count,t$effective_species.count)

pdf(paste0(outputPath,'raw counts vs effective species'))
plot(t$raw.count,t$effective_species.count, pch= c(22,21) ,col= c("10","4") ,main="raw counts vs effective species",xlab="raw counts",ylab="effective species count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
linewga <- lm(T.wga$effective_species.count~ T.wga$raw.count)
linenonwga <- lm(F.wga$effective_species.count~ F.wga$raw.count)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))
dev.off()

#just WGA 
#cor.test not working 
cor.test(T.wga$raw.count,T.wga$effective_species.count)
plot(T.wga$raw.count,T.wga$effective_species.count)
line<-lm(T.wga$effective_species.count~T.wga$raw.count)
abline(line)
anova(line)

#just non-WGA
cor.test(F.wga$raw.count,F.wga$effective_species.count)
plot(F.wga$raw.count,F.wga$effective_species.count)
line<-lm(F.wga$effective_species.count~F.wga$raw.count)
abline(line)
anova(line)

#raw read boxplot
t.test(T.wga$raw.count,F.wga$raw.count)

pdf(paste0(outputPath,'raw reads.pdf'))
boxplot(T.wga$raw.count, F.wga$raw.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="raw count",ylab="number of raw reads",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

boxplot(T.wga.ca$raw.count, T.wga.cf$raw.count,F.wga.ca$raw.count, F.wga.cf$raw.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" raw.count",ylab="read count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$raw.count,F.wga.ca$raw.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$raw.count,F.wga.cf$raw.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$raw.count,T.wga.cf$raw.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$raw.count,F.wga.cf$raw.count)

#bargraph 
ggplot(t,aes(pcrId,raw.count,fill=wga))+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#=========================================================================
#======================== usable reads - count ===========================
#=========================================================================

#usable counts vs effective species 
cor.test(t$usable.count,t$effective_species.count)
plot(t$usable.count,t$effective_species.count)

pdf(paste0(outputPath,'usable counts vs effective species'))
plot(t$usable.count,t$effective_species.count, pch= c(22,21) ,col= c("10","4") ,main="usable counts vs effective species",xlab="raw counts",ylab="effective species count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
linewga <- lm(T.wga$effective_species.count~ T.wga$usable.count)
linenonwga <- lm(F.wga$effective_species.count~ F.wga$usable.count)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))
dev.off()

#just WGA 
#cor.test not working 
cor.test(T.wga$usable.count,T.wga$effective_species.count)
plot(T.wga$usable.count,T.wga$effective_species.count)
line<-lm(T.wga$effective_species.count~T.wga$usable.count)
abline(line)
anova(line)

#just non-WGA
cor.test(F.wga$usable.count,F.wga$effective_species.count)
plot(F.wga$usable.count,F.wga$effective_species.count)
line<-lm(F.wga$effective_species.count~F.wga$usable.count)
abline(line)
anova(line)

#==== usable reads 
t.test(T.wga$usable.count,F.wga$usable.count)

pdf(paste0(outputPath,'usable seq.pdf'))
boxplot(T.wga$usable.count, F.wga$usable.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="usable seq",ylab="usable reads (%)",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()
#box plot for all fractions 
boxplot(T.wga.ca$usable.count, T.wga.cf$usable.count,F.wga.ca$usable.count, F.wga.cf$usable.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" usable.percent",ylab="percent",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$usable.count,F.wga.ca$usable.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$usable.count,F.wga.cf$usable.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$usable.count,T.wga.cf$usable.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$usable.count,F.wga.cf$usable.count)

#looking at correlations with usable reads (just non WGA samples )
#DNA yields 
cor.test(F.wga$usable.count,F.wga$dna)
plot(F.wga$usable.count,F.wga$dna)
line<-lm(F.wga$dna~F.wga$usable.count)
lineca<-lm(F.wga.ca$dna~F.wga.ca$usable.count)
abline(line)
abline(lineca)
anova(lineca)
anova(line)
#=====lymphocyte count 
cor.test(F.wga$usable.count,F.wga$lym.count)
plot(F.wga$usable.count,F.wga$lym.count)
line<-lm(F.wga$lym.count~F.wga$usable.count)
abline(line)
anova(line)

#====lymphocyte percentage
cor.test(F.wga$usable.count,F.wga$lym.percent)
plot(F.wga$usable.count,F.wga$lym.percent)
line<-lm(F.wga$lym.percent~F.wga$usable.count)
abline(line)
anova(line)

#bargraph 
ggplot(t,aes(pcrId,usable.count,fill=wga))+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))


#=========================================================================
#========================= usable reads - percent ========================
#=========================================================================

ggplot(t,aes(pcrId,usable.percent,fill=wga))+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))


#percent of usable reads and effective species 
cor.test(t$usable.percent,t$effective_species.count)
plot(t$usable.percent,t$effective_species.count)

pdf(paste0(outputPath,'usable seq percent vs effective species'))
plot(t$usable.percent,t$effective_species.count, pch= c(22,21) ,col= c("10","4") ,main="usable seq vs effective species",xlab="percent of usable seq",ylab="effective species count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
linewga <- lm(T.wga$effective_species.count~ T.wga$usable.percent)
linenonwga <- lm(F.wga$effective_species.count~ F.wga$usable.percent)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))
dev.off()

#just WGA 
#cor.test not working 
cor.test(T.wga$usable.percent,T.wga$effective_species.count)
plot(T.wga$usable.percent,T.wga$effective_species.count)
line<-lm(T.wga$effective_species.count~T.wga$usable.percent)
abline(line)
anova(line)

#just non-WGA
cor.test(F.wga$usable.percent,F.wga$effective_species.count)
plot(F.wga$usable.percent,F.wga$effective_species.count)
line<-lm(F.wga$effective_species.count~F.wga$usable.percent)
abline(line)
anova(line)

#== percent usable seq and clone counts 
cor.test(t$usable.percent,t$cln.count)
plot(t$usable.percent,t$cln.count)

pdf(paste0(outputPath,'usable seq percent vs number of unique clones'))
plot(t$usable.percent,t$cln.count, pch= c(22,21) ,col= c("10","4") ,main="usable seq vs clone count",xlab="percent of usable seq",ylab="unique clone count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
line<-lm(t$cln.count~t$usable.percent)
linewga <- lm(T.wga$cln.count~ T.wga$usable.percent)
linenonwga <- lm(F.wga$cln.count~ F.wga$usable.percent)
abline(line)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))
anova(line)

#just WGA 
#cor.test not working 
cor.test(T.wga$usable.percent,T.wga$cln.count)
plot(T.wga$usable.percent,T.wga$cln.count)
line<-lm(T.wga$cln.count~T.wga$usable.percent)
abline(line)
anova(line)

#just non-WGA
cor.test(F.wga$usable.percent,F.wga$cln.count)
plot(F.wga$usable.percent,F.wga$cln.count)
line<-lm(F.wga$cln.count~F.wga$usable.percent)
abline(line)
anova(line)

#boxplots 
t.test(T.wga$usable.percent,F.wga$usable.percent)

pdf(paste0(outputPath,'usable seq.pdf'))
boxplot(T.wga$usable.percent, F.wga$usable.percent,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="usable seq",ylab="usable reads (%)",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()
#box plot for all fractions 
boxplot(T.wga.ca$usable.percent, T.wga.cf$usable.percent,F.wga.ca$usable.percent, F.wga.cf$usable.percent,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" usable.percent",ylab="percent",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$usable.percent,F.wga.ca$usable.percent)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$usable.percent,F.wga.cf$usable.percent)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$usable.percent,T.wga.cf$usable.percent)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$usable.percent,F.wga.cf$usable.percent)

#DNA yields 
cor.test(F.wga$usable.percent,F.wga$dna)
plot(F.wga$usable.percent,F.wga$dna)
line<-lm(F.wga$dna~F.wga$usable.percent)
lineca<-lm(F.wga.ca$dna~F.wga.ca$usable.percent)
abline(line)
abline(lineca)
anova(lineca)
anova(line)
#=====lymphocyte count 
cor.test(F.wga$usable.percent,F.wga$lym.count)
plot(F.wga$usable.percent,F.wga$lym.count)
line<-lm(F.wga$lym.count~F.wga$usable.percent)
abline(line)
anova(line)

#====lymphocyte percentage
cor.test(F.wga$usable.percent,F.wga$lym.percent)
plot(F.wga$usable.percent,F.wga$lym.percent)
line<-lm(F.wga$lym.percent~F.wga$usable.percent)
abline(line)
anova(line)

#=========================================================================
#========================= uniq usable seq.count== ========================
#=========================================================================
#=====uniq usable seq.count
t.test(T.wga$uniq_usable_seq.count,F.wga$uniq_usable_seq.count)

pdf(paste0(outputPath,'usable seq.pdf'))
boxplot(T.wga$uniq_usable_seq.count, F.wga$uniq_usable_seq.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="usable seq",ylab="usable seq count",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()
#box plot for all fractions 
boxplot(T.wga.ca$uniq_usable_seq.count, T.wga.cf$uniq_usable_seq.count,F.wga.ca$uniq_usable_seq.count, F.wga.cf$uniq_usable_seq.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" usable.seq count",ylab="percent",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$uniq_usable_seq.count,F.wga.ca$uniq_usable_seq.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$uniq_usable_seq.count,F.wga.cf$uniq_usable_seq.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$uniq_usable_seq.count,T.wga.cf$uniq_usable_seq.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$uniq_usable_seq.count,F.wga.cf$uniq_usable_seq.count)

#looking at correlations with usable reads (just non WGA samples )
#===== DNA yields 
cor.test(F.wga$uniq_usable_seq.count,F.wga$dna)
plot(F.wga$uniq_usable_seq.count,F.wga$dna)
line<-lm(F.wga$dna~F.wga$uniq_usable_seq.count)
lineca<-lm(F.wga.ca$dna~F.wga.ca$uniq_usable_seq.count)
abline(line)
abline(lineca)
anova(lineca)
anova(line)
#=====lymphocyte count 
cor.test(F.wga$uniq_usable_seq.count,F.wga$lym.count)
plot(F.wga$uniq_usable_seq.count,F.wga$lym.count)
line<-lm(F.wga$lym.count~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#====lymphocyte percentage
cor.test(F.wga$uniq_usable_seq.count,F.wga$lym.percent)
plot(F.wga$uniq_usable_seq.count,F.wga$lym.percent)
line<-lm(F.wga$lym.percent~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)
#=========================================================================
#========================= cln- count== ========================
#=========================================================================
#====clonotype count 
t.test(T.wga$cln.count,F.wga$cln.count)

boxplot(T.wga$cln.count, F.wga$cln.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="clonotype count",ylab="cln.count",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))

#box plot for all fractions 
boxplot(T.wga.ca$cln.count, T.wga.cf$cln.count,F.wga.ca$cln.count, F.wga.cf$cln.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" clonotype count",ylab="cln.count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$cln.count,F.wga.ca$cln.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$cln.count,F.wga.cf$cln.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$cln.count,T.wga.cf$cln.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$cln.count,F.wga.cf$cln.count)

#=========================================================================
#========================= most popular cln- count== ========================
#=========================================================================
ggplot(t,aes(pcrId,w_most_pop_clns.count,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#=========================================================================
#========================= most popular cln- count percent== ========================
#=========================================================================
ggplot(t,aes(pcrId,w_most_pop_clns.percent,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

t.test(T.wga$w_most_pop_clns.percent,F.wga$w_most_pop_clns.percent)

boxplot(T.wga$w_most_pop_clns.percent, F.wga$w_most_pop_clns.percent,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="percent of the most popular clonotype",ylab="percent of total seq",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))

#box plot for all fractions 
boxplot(T.wga.ca$w_most_pop_clns.percent, T.wga.cf$w_most_pop_clns.percent,F.wga.ca$w_most_pop_clns.percent, F.wga.cf$w_most_pop_clns.percent,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" percent of the most popular clonotype",ylab="percent of total seq",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$w_most_pop_clns.percent,F.wga.ca$w_most_pop_clns.percent)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$w_most_pop_clns.percent,F.wga.cf$w_most_pop_clns.percent)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$w_most_pop_clns.percent,T.wga.cf$w_most_pop_clns.percent)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$w_most_pop_clns.percent,F.wga.cf$w_most_pop_clns.percent)

#=========================================================================
#==================== effective species - count ==========================
#=========================================================================

#t-test to compare effctive species overall
t.test (T.wga$effective_species.count,F.wga$effective_species.count)

pdf(paste0(outputPath,'boxplot effective species.pdf'))
boxplot(T.wga$effective_species.count, F.wga$effective_species.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="effective species",ylab="effective species",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()


#boxplot inculding all sample fractions 
pdf(paste0(outputPath,'boxplot effective species (comparing fractions).pdf'))
boxplot(T.wga.ca$effective_species.count, T.wga.cf$effective_species.count,F.wga.ca$effective_species.count, F.wga.cf$effective_species.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" effective species",ylab="# of effective species",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

  #T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$effective_species.count,F.wga.ca$effective_species.count)
  #T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$effective_species.count,F.wga.cf$effective_species.count)
  #T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$effective_species.count,T.wga.cf$effective_species.count)
  #T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$effective_species.count,F.wga.cf$effective_species.count)


#=========================================================================
#=================== alpha/shannon diversity - count =====================
#=========================================================================

t.test(T.wga$alpha_diversity.count,F.wga$alpha_diversity.count)

pdf(paste0(outputPath,'boxplot alpha diversity.pdf'))
boxplot(T.wga$alpha_diversity.count, F.wga$alpha_diversity.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="diversity",ylab="alpha diversity",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

#boxplot inculding all sample fractions 
boxplot(T.wga.ca$alpha_diversity.count, T.wga.cf$alpha_diversity.count,F.wga.ca$alpha_diversity.count, F.wga.cf$alpha_diversity.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" alpha diversity",ylab="alpha_diversity.count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))


#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$alpha_diversity.count,F.wga.ca$alpha_diversity.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$alpha_diversity.count,F.wga.cf$alpha_diversity.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$alpha_diversity.count,T.wga.cf$alpha_diversity.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$alpha_diversity.count,F.wga.cf$alpha_diversity.count)

#=========================================================================
#=================== diversity_normalized - count ========================
#=========================================================================

t.test(T.wga$diversity_normed.count,F.wga$diversity_normed.count)

boxplot(T.wga$diversity_normed.count, F.wga$diversity_normed.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="diversity",ylab="diversity_normed.count",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))

boxplot(T.wga.ca$diversity_normed.count, T.wga.cf$diversity_normed.count,F.wga.ca$diversity_normed.count, F.wga.cf$diversity_normed.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" diversity_normed.count",ylab="diversity_normed.count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$diversity_normed.count,F.wga.ca$diversity_normed.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$diversity_normed.count,F.wga.cf$diversity_normed.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$diversity_normed.count,T.wga.cf$diversity_normed.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$diversity_normed.count,F.wga.cf$diversity_normed.count)


#=========================================================================
#================== most popular clonotype percent =======================
#=========================================================================

t.test(T.wga$w_most_pop_clns.percent,F.wga$w_most_pop_clns.percent)

boxplot(T.wga$w_most_pop_clns.percent, F.wga$w_most_pop_clns.percent,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="percent of the most popular clonotype",ylab="percent of total seq",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))

#box plot for all fractions 
boxplot(T.wga.ca$w_most_pop_clns.percent, T.wga.cf$w_most_pop_clns.percent,F.wga.ca$w_most_pop_clns.percent, F.wga.cf$w_most_pop_clns.percent,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" percent of the most popular clonotype",ylab="percent of total seq",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$w_most_pop_clns.percent,F.wga.ca$w_most_pop_clns.percent)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$w_most_pop_clns.percent,F.wga.cf$w_most_pop_clns.percent)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$w_most_pop_clns.percent,T.wga.cf$w_most_pop_clns.percent)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$w_most_pop_clns.percent,F.wga.cf$w_most_pop_clns.percent)



