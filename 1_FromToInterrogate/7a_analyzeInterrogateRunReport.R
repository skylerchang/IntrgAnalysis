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
# 1 = csf run 1 = seq run 21
# 2 = csf run 2 = seq run 24 
# 3 = csf run 3 = seq run 26 
# 4 = excel file with both csf run 2 and 3 combined 
# 5 = excel file with all three seq runs 
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
  t$hemo.cells.L<-as.numeric(t$hemo.cells.L)
  t$hemo.cells.ul<-as.numeric(t$hemo.cells.ul)
  t$cytospin.count<-as.numeric(t$cytospin.count)
  t$lym.cyto.con<-as.numeric(t$lym.cyto.con)
  t$lym.cytocount<-as.numeric(t$lym.cytocount)
  t$lym.cytopercent<-as.numeric(t$lym.cytopercent)
  t$lym.hemo_con<-as.numeric(t$lym.hemo_con)
  t$hemolym.cells.ul<-as.numeric(t$hemolym.cells.ul)
  t$lym.hemo_count<-as.numeric(t$lym.hemo_count)
  t$pred.lym.cyto<-as.numeric(t$pred.lym.cyto)
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

#=============== for pdf output 
targetDir2<-'../Results/InterrogateRunReport/WGAgraphs/'

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

#remove run 1 for wga analysis (only do if analysisng WGA)
F.wga<- F.wga[-c(1:4,13:16,21:24,29:32,41:44,49:52),]

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


plot(t$raw.count,t$effective_species.count, pch= c(22,21) ,col= c("10","4") ,main="raw counts vs effective species",xlab="raw counts",ylab="effective species count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
linewga <- lm(T.wga$effective_species.count~ T.wga$raw.count)
linenonwga <- lm(F.wga$effective_species.count~ F.wga$raw.count)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))


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

#wga vs non-wga overall
pdf(paste0(targetDir2,'raw reads.pdf'))
boxplot(T.wga$raw.count, F.wga$raw.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="raw count",ylab="number of raw reads",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

#all samples 
pdf(paste0(targetDir2,'raw reads count (all samples).pdf'))
ggplot(t,aes(sampleId,raw.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'raw reads count (all samples no wga).pdf'))
ggplot(t,aes(sampleId,raw.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'raw reads count (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,raw.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

  #bargraph 
pdf(paste0(targetDir2,'raw read count (all samples- bargraph).pdf'))
ggplot(t,aes(wga,raw.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir2,'raw reads count (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,raw.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'raw read overall (non-wga).pdf'))
boxplot(F.wga.ca$raw.count, F.wga.cf$raw.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="Raw read count",ylab="Number of raw reads",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'raw reads count (wga).pdf'))
ggplot(T.wga,aes(sampleId,raw.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'raw read overall (wga).pdf'))
boxplot(T.wga.ca$raw.count, T.wga.cf$raw.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="Raw read count",ylab="Number of raw reads",xlab="Method",varwidth=TRUE)
dev.off()

#comparing all fractions for both wga and non-wga
pdf(paste0(targetDir2,'raw reads (all fractions).pdf'))
boxplot(T.wga.ca$raw.count, T.wga.cf$raw.count,F.wga.ca$raw.count, F.wga.cf$raw.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" raw.count",ylab="read count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

pdf(paste0(targetDir2,'raw reads (all fractions).pdf'))
boxplot(T.wga$raw.count,F.wga$raw.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "NON-WGA"),  cex.axis=1,main=" Raw read count",ylab="Number of raw reads",xlab="Method",varwidth=TRUE)
dev.off()

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


plot(t$usable.count,t$effective_species.count, pch= c(22,21) ,col= c("10","4") ,main="usable counts vs effective species",xlab="raw counts",ylab="effective species count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
linewga <- lm(T.wga$effective_species.count~ T.wga$usable.count)
linenonwga <- lm(F.wga$effective_species.count~ F.wga$usable.count)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))


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

#wga vs nonwga overall
pdf(paste0(targetDir2,'usable read count.pdf'))
boxplot(T.wga$usable.count, F.wga$usable.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="usable seq",ylab="usable reads (%)",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'usable read count (all samples).pdf'))
ggplot(t,aes(sampleId,usable.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'usable read count(all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,usable.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'usable read count (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,usable.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'usable read overall (non-wga).pdf'))
boxplot(F.wga.ca$usable.count, F.wga.cf$usable.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="Usable read count",ylab="Number of usable reads",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'usable read count (wga).pdf'))
ggplot(T.wga,aes(sampleId,usable.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'usable read overall (wga).pdf'))
boxplot(T.wga.ca$usable.count, T.wga.cf$usable.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="usable read count",ylab="number of usable reads",xlab="Method",varwidth=TRUE)
dev.off()

#box plot for all fractions 
pdf(paste0(targetDir2,'usable seq (all fractions).pdf'))
boxplot(T.wga.ca$usable.count, T.wga.cf$usable.count,F.wga.ca$usable.count, F.wga.cf$usable.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" usable.percent",ylab="percent",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

pdf(paste0(targetDir2,'raw reads (all fractions).pdf'))
boxplot(T.wga$usable.count,F.wga$usable.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "NON-WGA"),  cex.axis=1,main=" Usable read count",ylab="Number of usable reads",xlab="Method",varwidth=TRUE)
dev.off()

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
cor.test(F.wga$usable.count,F.wga$pre.dna)
plot(F.wga$usable.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$usable.count)
lineca<-lm(F.wga.ca$pre.dna~F.wga.ca$usable.count)
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
pdf(paste0(targetDir2,'usable read count (all samples- bargraph).pdf'))
ggplot(t,aes(wga,usable.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#=========================================================================
#========================= usable reads - percent ========================
#=========================================================================

ggplot(t,aes(pcrId,usable.percent,fill=wga))+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))


#percent of usable reads and effective species 
cor.test(t$usable.percent,t$effective_species.count)
plot(t$usable.percent,t$effective_species.count)


plot(t$usable.percent,t$effective_species.count, pch= c(22,21) ,col= c("10","4") ,main="usable seq vs effective species",xlab="percent of usable seq",ylab="effective species count")
legend ("topleft",inset=0.02,cex=0.7,title="Sample type",pch =c(22,21),col=c("10","4"),c("WGA","non-WGA"))
linewga <- lm(T.wga$effective_species.count~ T.wga$usable.percent)
linenonwga <- lm(F.wga$effective_species.count~ F.wga$usable.percent)
abline(linewga, col=c(10))
abline(linenonwga, col=c(4))

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

#overall wga vs nonwga
pdf(paste0(targetDir2,'usable seq percent.pdf'))
boxplot(T.wga$usable.percent, F.wga$usable.percent,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="usable seq",ylab="usable reads (%)",xlab="Method",varwidth=TRUE)
#legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'usable seq percent (all samples).pdf'))
ggplot(t,aes(sampleId,usable.percent,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'usable seq percent(all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,usable.percent,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

  #bargraph 
pdf(paste0(targetDir2,'usable seq percent (all samples- bargraph).pdf'))
ggplot(t,aes(wga,usable.percent,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir2,'usable seq percent (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,usable.percent,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'usable read overall percent(non-wga).pdf'))
boxplot(F.wga.ca$usable.percent, F.wga.cf$usable.percent,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="usable reads",ylab="percent of usable reads",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'usable seq percent (wga).pdf'))
ggplot(T.wga,aes(sampleId,usable.percent,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'usable read overall percent(wga).pdf'))
boxplot(T.wga.ca$usable.percent, T.wga.cf$usable.percent,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="usable reads",ylab="percent of usable reads",xlab="Method",varwidth=TRUE)
dev.off()

#box plot for all fractions 
pdf(paste0(targetDir2,'usable seq percent(all fractions).pdf'))
boxplot(T.wga.ca$usable.percent, T.wga.cf$usable.percent,F.wga.ca$usable.percent, F.wga.cf$usable.percent,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" usable.percent",ylab="percent",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

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
cor.test(F.wga$usable.percent,F.wga$pred.lym.cyto)
plot(F.wga$usable.percent,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$usable.percent)
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

pdf(paste0(targetDir2,'unique usable seq count.pdf'))
boxplot(T.wga$uniq_usable_seq.count, F.wga$uniq_usable_seq.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="usable seq",ylab="usable seq count",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'unique usable seq count.pdf (all samples).pdf'))
ggplot(t,aes(sampleId,uniq_usable_seq.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'unique usable seq count.pdf(all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,uniq_usable_seq.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

  #bargraph 
pdf(paste0(targetDir2,'unique usable seq count (all samples- bargraph).pdf'))
ggplot(t,aes(wga,uniq_usable_seq.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'unique usable seq count (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,uniq_usable_seq.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'unique usable read overall percent(non-wga).pdf'))
boxplot(F.wga.ca$uniq_usable_seq.count, F.wga.cf$uniq_usable_seq.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="unique usable sequences",ylab="seq count",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'unique usable seq count (wga).pdf'))
ggplot(T.wga,aes(sampleId,uniq_usable_seq.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'unique usable read overall percent(wga).pdf'))
boxplot(T.wga.ca$uniq_usable_seq.count, T.wga.cf$uniq_usable_seq.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="unique usable sequences",ylab="seq count",xlab="Method",varwidth=TRUE)
dev.off()

#box plot for all fractions 
pdf(paste0(targetDir2,'unique usable seq count(all fractions).pdf'))
boxplot(T.wga.ca$uniq_usable_seq.count, T.wga.cf$uniq_usable_seq.count,F.wga.ca$uniq_usable_seq.count, F.wga.cf$uniq_usable_seq.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" usable.seq count",ylab="count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

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

# bar graph 
ggplot(t,aes(pcrId,uniq_usable_seq.count,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#=========================================================================
#========================= cln- count== ========================
#=========================================================================
#====clonotype count 
t.test(T.wga$cln.count,F.wga$cln.count)

pdf(paste0(targetDir2,'clonotype count.pdf'))
boxplot(T.wga$cln.count, F.wga$cln.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="clonotype count",ylab="cln.count",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'clonotype count (all samples).pdf'))
ggplot(t,aes(sampleId,cln.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'clonotype count (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,cln.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+xlab("WGA method")+ ylab("Clonotype count") + ggtitle("Clonotype count") +theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# using dog id instead of submission number 
pdf(paste0(targetDir2,'clonotype count (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,cln.count,fill=sampleFraction))+facet_grid(.~dog.id)+geom_boxplot()+xlab("WGA method")+ ylab("Clonotype count") + ggtitle("Clonotype count") +scale_fill_discrete( name = "Sample Fraction",labels = c("ca", "cf"))+ scale_x_discrete( name = "WGA method",labels = c("non-WGA", "WGA"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


pdf(paste0(targetDir2,'clonotype count (all samples by wga and fraction).pdf'))
ggplot(t,aes(dog.id,cln.count,fill=sampleFraction))+facet_grid(.~wga)+geom_boxplot()+xlab("WGA method")+ ylab("Clonotype count") + ggtitle("Clonotype count") +scale_fill_discrete( name = "Sample Fraction",labels = c("ca", "cf"))+ scale_x_discrete( name = "Patient")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

  #bargraph 
pdf(paste0(targetDir2,'clonotype count (all samples- bargraph).pdf'))
ggplot(t,aes(wga,cln.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity") + scale_x_discrete( name = "Sample Fraction",labels = c("ca", "cf"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'cln clonotype count (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,cln.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'clonotype count overall (non-wga).pdf'))
boxplot(F.wga.ca$cln.count, F.wga.cf$cln.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="clonotype",ylab="cln count",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'cln clonotype count (wga).pdf'))
ggplot(T.wga,aes(sampleId,cln.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'clonotype count overall (wga).pdf'))
boxplot(T.wga.ca$cln.count, T.wga.cf$cln.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="clonotype",ylab="cln count",xlab="Method",varwidth=TRUE)
dev.off()

#box plot for all fractions 
pdf(paste0(targetDir2,'clonotype count (all fractions).pdf'))
boxplot(T.wga.ca$cln.count, T.wga.cf$cln.count,F.wga.ca$cln.count, F.wga.cf$cln.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" clonotype count",ylab="cln.count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$cln.count,F.wga.ca$cln.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$cln.count,F.wga.cf$cln.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$cln.count,T.wga.cf$cln.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$cln.count,F.wga.cf$cln.count)

# bar graph 
ggplot(t,aes(pcrId,cln.count,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#=========================================================================
#========================= most popular cln- count== ========================
#=========================================================================
ggplot(t,aes(pcrId,w_most_pop_clns.count,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

t.test(T.wga$w_most_pop_clns.count,F.wga$w_most_pop_clns.percent)

#boxplot overall wga vs nonwga
pdf(paste0(targetDir2,'most abundant cln clonotype count.pdf'))
boxplot(T.wga$w_most_pop_clns.count, F.wga$w_most_pop_clns.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="percent of the most popular clonotype",ylab="percent of total seq",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'most abundant cln clonotype count(all samples).pdf'))
ggplot(t,aes(sampleId,w_most_pop_clns.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off() 
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'most abundant cln clonotype count (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,w_most_pop_clns.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

  #bargraph 
pdf(paste0(targetDir2,'most abundant cln clonotype count (all samples- bargraph).pdf'))
ggplot(t,aes(wga,w_most_pop_clns.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'most abundant cln clonotype count (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,w_most_pop_clns.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'most abundant cln clonotype count overall (non-wga).pdf'))
boxplot(F.wga.ca$w_most_pop_clns.count, F.wga.cf$w_most_pop_clns.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="count of the most popular clonotype",ylab="cln count",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'most abundant cln clonotype count (wga).pdf'))
ggplot(T.wga,aes(sampleId,w_most_pop_clns.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'most abundant cln clonotype count overall (wga).pdf'))
boxplot(T.wga.ca$w_most_pop_clns.count, T.wga.cf$w_most_pop_clns.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="count of the most popular clonotype",ylab="cln count",xlab="Method",varwidth=TRUE)
dev.off()

#box plot for all fractions 
pdf(paste0(targetDir2,'most abundant cln clonotype count (all fractions).pdf'))
boxplot(T.wga.ca$w_most_pop_clns.count, T.wga.cf$w_most_pop_clns.count,F.wga.ca$w_most_pop_clns.percent, F.wga.cf$w_most_pop_clns.percent,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" percent of the most popular clonotype",ylab="percent of total seq",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$w_most_pop_clns.count,F.wga.ca$w_most_pop_clns.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$w_most_pop_clns.count,F.wga.cf$w_most_pop_clns.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$w_most_pop_clns.count,T.wga.cf$w_most_pop_clns.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$w_most_pop_clns.count,F.wga.cf$w_most_pop_clns.count)

#=========================================================================
#========================= most popular cln- count percent== ========================
#=========================================================================
ggplot(t,aes(pcrId,w_most_pop_clns.percent,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

t.test(T.wga$w_most_pop_clns.percent,F.wga$w_most_pop_clns.percent)

#boxplot overall wga vs nonwga
pdf(paste0(targetDir2,'most abundant cln clonotype percent.pdf'))
boxplot(T.wga$w_most_pop_clns.percent, F.wga$w_most_pop_clns.percent,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="percent of the most popular clonotype",ylab="percent of total seq",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'most abundant cln clonotype percent(all samples).pdf'))
ggplot(t,aes(sampleId,w_most_pop_clns.percent,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off() 
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'most abundant cln clonotype percent (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,w_most_pop_clns.percent,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

  #bargraph 
pdf(paste0(targetDir2,'most abundant cln clonotype percent (all samples- bargraph).pdf'))
ggplot(t,aes(wga,w_most_pop_clns.percent,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'most abundant cln clonotype percent (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,w_most_pop_clns.percent,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'most abundant cln clonotype percent overall (non-wga).pdf'))
boxplot(F.wga.ca$w_most_pop_clns.percent, F.wga.cf$w_most_pop_clns.percent,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="portion of the most popular clonotype",ylab="percent",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'most abundant cln clonotype percent (wga).pdf'))
ggplot(T.wga,aes(sampleId,w_most_pop_clns.percent,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'most abundant cln clonotype percent overall (wga).pdf'))
boxplot(T.wga.ca$w_most_pop_clns.percent, T.wga.cf$w_most_pop_clns.percent,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="portion of the most popular clonotype",ylab="percent",xlab="Method",varwidth=TRUE)
dev.off()

#box plot for all fractions 
pdf(paste0(targetDir2,'most abundant cln clonotype percent (all fractions).pdf'))
boxplot(T.wga.ca$w_most_pop_clns.percent, T.wga.cf$w_most_pop_clns.percent,F.wga.ca$w_most_pop_clns.percent, F.wga.cf$w_most_pop_clns.percent,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" percent of the most popular clonotype",ylab="percent of total seq",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

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

pdf(paste0(targetDir2,'boxplot effective species.pdf'))
boxplot(T.wga$effective_species.count, F.wga$effective_species.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="effective species",ylab="effective species",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'effective species (all samples).pdf'))
ggplot(t,aes(sampleId,effective_species.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off() 
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'effective species (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,effective_species.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+xlab("WGA method")+ ylab("Number of effective species") + ggtitle("Effective species count") +theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
# using dog id instead of submission number 
pdf(paste0(targetDir2,'effective species (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,effective_species.count,fill=sampleFraction))+facet_grid(.~dog.id)+geom_boxplot()+xlab("WGA method")+ ylab("Effective species count") + ggtitle("Effective species count") +scale_fill_discrete( name = "Sample Fraction",labels = c("ca", "cf"))+ scale_x_discrete( name = "WGA method",labels = c("non-WGA", "WGA"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
  #bargraph 
pdf(paste0(targetDir2,'effective species (all samples- bargraph).pdf'))
ggplot(t,aes(wga,effective_species.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'boxplot effective species (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,effective_species.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'effective species overall (non-wga).pdf'))
boxplot(F.wga.ca$effective_species.count, F.wga.cf$effective_species.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="effective species",ylab="count",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'boxplot effective species (wga).pdf'))
ggplot(T.wga,aes(sampleId,effective_species.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'effective species overall (wga).pdf'))
boxplot(T.wga.ca$effective_species.count, T.wga.cf$effective_species.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="effective species",ylab="count",xlab="Method",varwidth=TRUE)
dev.off()

#boxplot inculding all sample fractions 
pdf(paste0(targetDir2,'boxplot effective species (comparing fractions).pdf'))
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

pdf(paste0(targetDir2,'boxplot alpha diversity.pdf'))
boxplot(T.wga$alpha_diversity.count, F.wga$alpha_diversity.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="diversity",ylab="alpha diversity",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'alpha diversity (all samples).pdf'))
ggplot(t,aes(sampleId,alpha_diversity.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'alpha diversity (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,alpha_diversity.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#bargraph 
pdf(paste0(targetDir2,'alpha diversity  (all samples- bargraph).pdf'))
ggplot(t,aes(wga,alpha_diversity.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'boxplot alpha diversity (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,alpha_diversity.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'alpha diveraity overall (non-wga).pdf'))
boxplot(F.wga.ca$alpha_diversity.count, F.wga.cf$alpha_diversity.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="alpha diversity",ylab="count",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'boxplot alpha diversity (wga).pdf'))
ggplot(T.wga,aes(sampleId,alpha_diversity.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'alpha diveraity overall (wga).pdf'))
boxplot(T.wga.ca$alpha_diversity.count, T.wga.cf$alpha_diversity.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="alpha diversity",ylab="count",xlab="Method",varwidth=TRUE)
dev.off()

#boxplot inculding all sample fractions 
pdf(paste0(targetDir2,'boxplot alpha diversity (all fractions).pdf'))
boxplot(T.wga.ca$alpha_diversity.count, T.wga.cf$alpha_diversity.count,F.wga.ca$alpha_diversity.count, F.wga.cf$alpha_diversity.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" alpha diversity",ylab="alpha_diversity.count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

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

#overall
pdf(paste0(targetDir2,'normalized diversity.pdf'))
boxplot(T.wga$diversity_normed.count, F.wga$diversity_normed.count,col=c("10","4"),outcol=c("10","4"), names=c("WGA", "non-WGA"),  cex.axis=1,main="diversity",ylab="diversity_normed.count",xlab="Method",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","non-WGA"))
dev.off()

pdf(paste0(targetDir2,'normalized diversity (all samples).pdf'))
ggplot(t,aes(sampleId,diversity_normed.count,fill=wga))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
  #boxplot with more then one layer of organization 
pdf(paste0(targetDir2,'normalized diversity (all samples by wga and fraction).pdf'))
ggplot(t,aes(wga,diversity_normed.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#bargraph 
pdf(paste0(targetDir2,'normalized diversity (all samples- bargraph).pdf'))
ggplot(t,aes(wga,diversity_normed.count,fill=sampleFraction))+facet_grid(.~submissionId)+geom_bar(position = "dodge", stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#just non-wga 
pdf(paste0(targetDir,'normalized diversity (non-wga).pdf'))
ggplot(F.wga,aes(sampleId,diversity_normed.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir,'normalized diveraity overall (non-wga).pdf'))
boxplot(F.wga.ca$diversity_normed.count, F.wga.cf$diversity_normed.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="normalized diversity",ylab="count",xlab="Method",varwidth=TRUE)
dev.off()

#just WGA
pdf(paste0(targetDir2,'normalized diversity (wga).pdf'))
ggplot(T.wga,aes(sampleId,diversity_normed.count,fill=sampleFraction))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf(paste0(targetDir2,'normalized diveraity overall (non-wga).pdf'))
boxplot(T.wga.ca$diversity_normed.count, T.wga.cf$diversity_normed.count,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="normalized diversity",ylab="count",xlab="Method",varwidth=TRUE)
dev.off()

#comparing all fractions 
pdf(paste0(targetDir2,'normalized diversity (all fractions).pdf'))
boxplot(T.wga.ca$diversity_normed.count, T.wga.cf$diversity_normed.count,F.wga.ca$diversity_normed.count, F.wga.cf$diversity_normed.count,col=c("10","10","4","4"),outcol=c("10","10","4","4"), names=c("Cell-associated", "Cell-free","Cell-associated", "Cell-free"),  cex.axis=1,main=" diversity_normed.count",ylab="diversity_normed.count",xlab="Sample fraction",varwidth=TRUE)
legend ("topleft",inset=0.02,cex=1,title="Sample type",pch =c(21,21),col=c("10","4"),c("WGA","Non-WGA"))
dev.off()

#T-TEST wga ca fraction vs non-wga cf fraction 
t.test (T.wga.ca$diversity_normed.count,F.wga.ca$diversity_normed.count)
#T-TEST wga cf fraction vs non-wga cf fraction 
t.test (T.wga.cf$diversity_normed.count,F.wga.cf$diversity_normed.count)
#T-TEST wga ca fraction vs wga cf fraction 
t.test (T.wga.ca$diversity_normed.count,T.wga.cf$diversity_normed.count)
#T-TEST non-wga ca fraction vs non-wga ca fraction 
t.test (F.wga.ca$diversity_normed.count,F.wga.cf$diversity_normed.count)

# bar graph 
ggplot(t,aes(pcrId,diversity_normed.count,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#looking at the confidence of the normalized diversity count 
# plot normalized diversity vs total read count 

pdf(paste0(targetDir,'normalized diversity vs total read count (by WGA treatment).pdf'))
ggplot(t,aes(raw.count,diversity_normed.count,fill=wga))+geom_point(size=2, shape=23)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# by seq runs 
pdf(paste0(targetDir,'normalized diversity vs total read count(by runs).pdf'))
ggplot(t,aes(raw.count,diversity_normed.count,fill=seq.run))+geom_point(size=2, shape=23)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#non-wga
pdf(paste0(targetDir,'normalized diversity vs total read count(nonwga).pdf'))
plot(F.wga$raw.count,F.wga$diversity_normed.count,main="normalized diversity vs total read count",xlab="raw counts",ylab="normalized diversity count")
dev.off()

# wga 
pdf(paste0(targetDir,'normalized diversity vs total read count (WGA).pdf'))
plot(T.wga$raw.count,T.wga$diversity_normed.count,main="normalized diversity vs total read count",xlab="raw counts",ylab="normalized diversity count")
dev.off()


#=================================================================================
#========================= determining perdicted variables =======================
#=================================================================================

#first we will create two new variables one for lym count and one for lym concentration
# this will use cyto values for samples below 200 cytospin count and hemo values for samples with a cytospin count of 200

#create 'final.lymcount' 
t$final.lymcount<-t$lym.cytocount
t[t$cytospin.count=="200" & !is.na(t$lym.cytocount),"final.lymcount"]<-t[t$cytospin.count == "200" & !is.na(t$lym.cytocount),"lym.hemo_count"]
t[t$cytospin.count=="100" & !is.na(t$lym.cytocount),"final.lymcount"]<-t[t$cytospin.count == "100" & !is.na(t$lym.cytocount),"lym.hemo_count"]

#create 'final. lym concentration (cells/ul) ->
t$final.lymcon<-t$lym.cyto_con
t[t$cytospin.count=="200" & !is.na(t$lym.cyto_con),"final.lymcon"]<-t[t$cytospin.count == "200" & !is.na(t$lym.cyto_con),"hemolym.cells.ul"]
t[t$cytospin.count=="100" & !is.na(t$lym.cyto_con),"final.lymcon"]<-t[t$cytospin.count == "100" & !is.na(t$lym.cyto_con),"hemolym.cells.ul"]

#========== analysis 
#========== deteriming if there any perdicted variables to determine if a sample is ideal for NGS
#variables to look at 
#1. DNA yields
#2. hemo (cell concentration(cells/ul))
#3. hemo (lym concentration (cells/ul))
#4. hemo (lym count (calculated for 200 ul same as the cytospin))
#5. cyto (lym concentration (cells/ul)) 
#6. cyto (lym count)
#7. using final lym concentration 
#8. using final lym count 
#9. predicted cyto spin count 

#==============1. DNA yields 
#============================= need to fix 
pdf(paste0(targetDir2,'DNA yileds (all samples).pdf'))
ggplot(t,aes(sampleId,postwga.dna,fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
dev.off()

pdf(paste0(targetDir2,'DNA yileds (input vs output).pdf'))
boxplot(F.wga$postwga.dna, T.wga$postwga.dna,col=c("10","4"),outcol=c("10","4"), names=c("input", "output"),  cex.axis=1,main="WGA DNA yields",ylab="DNA yields (ng/ul)",varwidth=TRUE)
dev.off()

pdf(paste0(targetDir,'DNA yileds (overall).pdf'))
boxplot(F.wga$pre.dna,col=c("10"),outcol=c("10"), names=c("All fractions"),  cex.axis=1,main="DNA yields",ylab="DNA yields (ng/ul)",varwidth=TRUE)
dev.off()
summary (F.wga$pre.dna)

pdf(paste0(targetDir2,'DNA yileds (ca vs cf).pdf'))
boxplot(F.wga.ca$pre.dna, F.wga.cf$pre.dna,col=c("10","4"),outcol=c("10","4"), names=c("cell-associated", "cell-free"),  cex.axis=1,main="DNA yields",ylab="DNA yields (ng/ul)",varwidth=TRUE)
dev.off()
summary (F.wga.ca$pre.dna)
summary (F.wga.cf$pre.dna)

t.test (F.wga.ca$pre.dna,F.wga.cf$pre.dna)
#bwtween runs 
pdf(paste0(targetDir,'DNA yields (btw runs).pdf'))
ggplot(F.wga,aes(sampleFraction,pre.dna,fill=sampleFraction))+facet_grid(.~seq.run)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#correlations 
#raw reads
cor.test(F.wga$raw.count,F.wga$pre.dna)
plot(F.wga$raw.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$raw.count)
abline(line)
anova(line)

cor.test(F.wga.ca$raw.count,F.wga.ca$pre.dna)
cor.test(F.wga.cf$raw.count,F.wga.cf$pre.dna)

#usable reads percent
cor.test(F.wga$usable.percent,F.wga$pre.dna)
plot(F.wga$usable.percent,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$usable.percent)
abline(line)
anova(line)


#usable reads count
cor.test(F.wga$usable.count,F.wga$pre.dna)
plot(F.wga$usable.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$usable.count)
abline(line)
anova(line)

cor.test(F.wga.ca$usable.count,F.wga.ca$pre.dna)
cor.test(F.wga.cf$usable.count,F.wga.cf$pre.dna)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$pre.dna)
plot(F.wga$uniq_usable_seq.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$pre.dna)
plot(F.wga$cln.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$cln.count)
abline(line)
anova(line)

cor.test(F.wga.ca$cln.count,F.wga.ca$pre.dna)
cor.test(F.wga.cf$cln.count,F.wga.cf$pre.dna)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$pre.dna)
plot(F.wga$alpha_diversity.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$alpha_diversity.count)
abline(line)
anova(line)

cor.test(F.wga.ca$alpha_diversity.count,F.wga.ca$pre.dna)
cor.test(F.wga.cf$alpha_diversity.count,F.wga.cf$pre.dna)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$pre.dna)
plot(F.wga$diversity_normed.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$diversity_normed.count)
abline(line)
anova(line)


#effective species 
cor.test(F.wga$effective_species.count,F.wga$pre.dna)
plot(F.wga$effective_species.count,F.wga$pre.dna)
line<-lm(F.wga$pre.dna~F.wga$effective_species.count)
abline(line)
anova(line)

cor.test(F.wga.ca$effective_species.count,F.wga.ca$pre.dna)
cor.test(F.wga.cf$effective_species.count,F.wga.cf$pre.dna)
#==============2. hemo (cell concentration(cells/ul))
#====================================================
#bargraph
ggplot(t,aes(submissionId,hemo.cells.ul))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$hemo.cells.ul)
plot(F.wga$raw.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$raw.count)
abline(line)
anova(line)

cor.test(F.wga.ca$raw.count,F.wga.ca$hemo.cells.ul)
cor.test(F.wga.cf$raw.count,F.wga.cf$hemo.cells.ul)

#usable reads percent
cor.test(F.wga$usable.percent,F.wga$hemo.cells.ul)
plot(F.wga$usable.percent,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$hemo.cells.ul)
plot(F.wga$usable.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$usable.count)
abline(line)
anova(line)

cor.test(F.wga.ca$usable.count,F.wga.ca$hemo.cells.ul)
cor.test(F.wga.cf$usable.count,F.wga.cf$hemo.cells.ul)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$hemo.cells.ul)
plot(F.wga$uniq_usable_seq.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$hemo.cells.ul)
plot(F.wga$cln.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$cln.count)
abline(line)
anova(line)

cor.test(F.wga.ca$cln.count,F.wga.ca$hemo.cells.ul)
cor.test(F.wga.cf$cln.count,F.wga.cf$hemo.cells.ul)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$hemo.cells.ul)
plot(F.wga$alpha_diversity.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$alpha_diversity.count)
abline(line)
anova(line)

cor.test(F.wga.ca$alpha_diversity.count,F.wga.ca$hemo.cells.ul)
cor.test(F.wga.cf$alpha_diversity.count,F.wga.cf$hemo.cells.ul)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$hemo.cells.ul)
plot(F.wga$diversity_normed.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$hemo.cells.ul)
plot(F.wga$effective_species.count,F.wga$hemo.cells.ul)
line<-lm(F.wga$hemo.cells.ul~F.wga$effective_species.count)
abline(line)
anova(line)

cor.test(F.wga.ca$effective_species.count,F.wga.ca$hemo.cells.ul)
cor.test(F.wga.cf$effective_species.count,F.wga.cf$hemo.cells.ul)
#==============3. hemo (lym concentration (cells/ul))
#====================================================
#bargraph
ggplot(t,aes(submissionId,hemolym.cells.ul))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$hemolym.cells.ul)
plot(F.wga$raw.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$raw.count)
abline(line)
anova(line)

#usable reads percent 
cor.test(F.wga$usable.percent,F.wga$hemolym.cells.ul)
plot(F.wga$usable.percent,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$hemolym.cells.ul)
plot(F.wga$usable.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$usable.count)
abline(line)
anova(line)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$hemolym.cells.ul)
plot(F.wga$uniq_usable_seq.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$hemolym.cells.ul)
plot(F.wga$cln.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$cln.count)
abline(line)
anova(line)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$hemolym.cells.ul)
plot(F.wga$alpha_diversity.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$alpha_diversity.count)
abline(line)
anova(line)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$hemolym.cells.ul)
plot(F.wga$diversity_normed.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$hemolym.cells.ul)
plot(F.wga$effective_species.count,F.wga$hemolym.cells.ul)
line<-lm(F.wga$hemolym.cells.ul~F.wga$effective_species.count)
abline(line)
anova(line)



#==============4. hemo (lym count (calculated for 200 ul same as the cytospin))
#==============================================================================
#bargraph
ggplot(t,aes(submissionId,lym.hemo_count))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 


#raw reads
cor.test(F.wga$raw.count,F.wga$lym.hemo_count)
plot(F.wga$raw.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$raw.count)
abline(line)
anova(line)

#usable reads percent 
cor.test(F.wga$usable.percent,F.wga$lym.hemo_count)
plot(F.wga$usable.percent,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$lym.hemo_count)
plot(F.wga$usable.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$usable.count)
abline(line)
anova(line)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$lym.hemo_count)
plot(F.wga$uniq_usable_seq.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$lym.hemo_count)
plot(F.wga$cln.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$cln.count)
abline(line)
anova(line)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$lym.hemo_count)
plot(F.wga$alpha_diversity.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$alpha_diversity.count)
abline(line)
anova(line)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$lym.hemo_count)
plot(F.wga$diversity_normed.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$lym.hemo_count)
plot(F.wga$effective_species.count,F.wga$lym.hemo_count)
line<-lm(F.wga$lym.hemo_count~F.wga$effective_species.count)
abline(line)
anova(line)


#==============5. cyto (lym concentration (cells/ul)) 
#====================================================
#bargraph
ggplot(t,aes(submissionId,lym.cyto_con))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$lym.cyto_con)
plot(F.wga$raw.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$raw.count)
abline(line)
anova(line)

#usable reads percent 
cor.test(F.wga$usable.percent,F.wga$lym.cyto_con)
plot(F.wga$usable.percent,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$lym.cyto_con)
plot(F.wga$usable.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$usable.count)
abline(line)
anova(line)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$lym.cyto_con)
plot(F.wga$uniq_usable_seq.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$lym.cyto_con)
plot(F.wga$cln.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$cln.count)
abline(line)
anova(line)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$lym.cyto_con)
plot(F.wga$alpha_diversity.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$alpha_diversity.count)
abline(line)
anova(line)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$lym.cyto_con)
plot(F.wga$diversity_normed.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$lym.cyto_con)
plot(F.wga$effective_species.count,F.wga$lym.cyto_con)
line<-lm(F.wga$lym.cyto_con~F.wga$effective_species.count)
abline(line)
anova(line)


#==============6. cyto (lym count) 
#=================================
#bargraph
ggplot(t,aes(submissionId,lym.cytocount))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$lym.cytocount)
plot(F.wga$raw.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$raw.count)
abline(line)
anova(line)

#usable reads percent 
cor.test(F.wga$usable.percent,F.wga$lym.cytocount)
plot(F.wga$usable.percent,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$lym.cytocount)
plot(F.wga$usable.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$usable.count)
abline(line)
anova(line)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$lym.cytocount)
plot(F.wga$uniq_usable_seq.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$lym.cytocount)
plot(F.wga$cln.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$cln.count)
abline(line)
anova(line)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$lym.cytocount)
plot(F.wga$alpha_diversity.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$alpha_diversity.count)
abline(line)
anova(line)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$lym.cytocount)
plot(F.wga$diversity_normed.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$lym.cytocount)
plot(F.wga$effective_species.count,F.wga$lym.cytocount)
line<-lm(F.wga$lym.cytocount~F.wga$effective_species.count)
abline(line)
anova(line)

#==============7. using final lym concentration 
#==============================================
# comparing concentration from cytospin vs hemocytometer 

boxplot(F.wga$lym.cyto_con, F.wga$hemolym.cells.ul,col=c("10","4"),outcol=c("10","4"), names=c("cytospin", "hemocytometer"),  cex.axis=1,main="lymphocyte concentration",ylab="cells/ul",xlab="cell count method",varwidth=TRUE)

t.test (F.wga$lym.cyto_con,F.wga$hemolym.cells.ul)

#bargraph
ggplot(t,aes(submissionId,final.lymcon))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$final.lymcon)
plot(F.wga$raw.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$raw.count)
abline(line)
anova(line)

#usable reads percent 
cor.test(F.wga$usable.percent,F.wga$final.lymcon)
plot(F.wga$usable.percent,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$final.lymcon)
plot(F.wga$usable.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$usable.count)
abline(line)
anova(line)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$final.lymcon)
plot(F.wga$uniq_usable_seq.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$final.lymcon)
plot(F.wga$cln.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$cln.count)
abline(line)
anova(line)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$final.lymcon)
plot(F.wga$alpha_diversity.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$alpha_diversity.count)
abline(line)
anova(line)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$final.lymcon)
plot(F.wga$diversity_normed.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$final.lymcon)
plot(F.wga$effective_species.count,F.wga$final.lymcon)
line<-lm(F.wga$final.lymcon~F.wga$effective_species.count)
abline(line)
anova(line)


#==============8. using final lym count 
#==============================================
# comparing lym counts from cytospin vs hemocytometer (for 200 ul)

boxplot(F.wga$lym.cytocount, F.wga$lym.hemo_count,col=c("10","4"),outcol=c("10","4"), names=c("cytospin", "hemocytometer"),  cex.axis=1,main="lymphocyte count",ylab="cells",xlab="cell count method",varwidth=TRUE)

t.test (F.wga$lym.cytocount,F.wga$lym.hemo_count)

#bargraph
ggplot(t,aes(submissionId,final.lymcount))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$final.lymcount)
plot(F.wga$raw.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$raw.count)
abline(line)
anova(line)

#usable reads percent 
cor.test(F.wga$usable.percent,F.wga$final.lymcount)
plot(F.wga$usable.percent,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$final.lymcount)
plot(F.wga$usable.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$usable.count)
abline(line)
anova(line)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$final.lymcount)
plot(F.wga$uniq_usable_seq.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$final.lymcount)
plot(F.wga$cln.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$cln.count)
abline(line)
anova(line)

#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$final.lymcount)
plot(F.wga$alpha_diversity.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$alpha_diversity.count)
abline(line)
anova(line)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$final.lymcount)
plot(F.wga$diversity_normed.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$final.lymcount)
plot(F.wga$effective_species.count,F.wga$final.lymcount)
line<-lm(F.wga$final.lymcount~F.wga$effective_species.count)
abline(line)
anova(line)

#=================================================================================
#predicted cyto
#=================================================================================
#bargraph
ggplot(t,aes(submissionId,pred.lym.cyto))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#correlations with seq variables 

#raw reads
cor.test(F.wga$raw.count,F.wga$pred.lym.cyto)
pdf(paste0(targetDir2,'Raw reads vs lymphocyte abundace .pdf'))
plot(F.wga$pred.lym.cyto,F.wga$raw.count,xlab ="Lymphocyte count",ylab= "Raw read count")
line<-lm(F.wga$raw.count~F.wga$pred.lym.cyto)
abline(line)
dev.off()
anova(line)

pdf(paste0(targetDir2,'Raw reads vs lymphocyte abundace .pdf'))
ggplot(F.wga,aes(pred.lym.cyto,raw.count,fill=sampleFraction))+geom_point(size=2, shape=23)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

cor.test(F.wga.ca$raw.count,F.wga.ca$pred.lym.cyto)
cor.test(F.wga.cf$raw.count,F.wga.cf$pred.lym.cyto)

#usable reads percent
cor.test(F.wga$usable.percent,F.wga$pred.lym.cyto)
plot(F.wga$usable.percent,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$usable.percent)
abline(line)
anova(line)

#usable reads count
cor.test(F.wga$usable.count,F.wga$pred.lym.cyto)
plot(F.wga$pred.lym.cyto,F.wga$usable.count,main= "Predicted variable",xlab ="Lymphocyte count",ylab= "Usable read count")
line<-lm(F.wga$usable.count~F.wga$pred.lym.cyto)
abline(line)
anova(line)

pdf(paste0(targetDir2,'usable reads vs lymphocyte abundace .pdf'))
ggplot(F.wga,aes(pred.lym.cyto,usable.count,fill=sampleFraction))+geom_point(size=2, shape=23)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

cor.test(F.wga.ca$usable.count,F.wga.ca$pred.lym.cyto)
cor.test(F.wga.cf$usable.count,F.wga.cf$pred.lym.cyto)

#usable seq 
cor.test(F.wga$uniq_usable_seq.count,F.wga$pred.lym.cyto)
plot(F.wga$uniq_usable_seq.count,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$uniq_usable_seq.count)
abline(line)
anova(line)

#clonotype count
cor.test(F.wga$cln.count,F.wga$pred.lym.cyto)
plot(F.wga$cln.count,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$cln.count)
abline(line)
anova(line)

cor.test(F.wga.ca$cln.count,F.wga.ca$pred.lym.cyto)
cor.test(F.wga.cf$cln.count,F.wga.cf$pred.lym.cyto)


#alpha diversity
cor.test(F.wga$alpha_diversity.count,F.wga$pred.lym.cyto)
plot(F.wga$alpha_diversity.count,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$alpha_diversity.count)
abline(line)
anova(line)

cor.test(F.wga.ca$alpha_diversity.count,F.wga.ca$pred.lym.cyto)
cor.test(F.wga.cf$alpha_diversity.count,F.wga.cf$pred.lym.cyto)

#normalized diversity
cor.test(F.wga$diversity_normed.count,F.wga$pred.lym.cyto)
plot(F.wga$diversity_normed.count,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$diversity_normed.count)
abline(line)
anova(line)

#effective species 
cor.test(F.wga$effective_species.count,F.wga$pred.lym.cyto)
plot(F.wga$effective_species.count,F.wga$pred.lym.cyto)
line<-lm(F.wga$pred.lym.cyto~F.wga$effective_species.count)
abline(line)
anova(line)

cor.test(F.wga.ca$effective_species.count,F.wga.ca$pred.lym.cyto)
cor.test(F.wga.cf$effective_species.count,F.wga.cf$pred.lym.cyto)
#=================================================================================
#=============================== sample Hierarchy ================================
#=================================================================================

# three diversity variables 
#1. alpha diversity
#2. normalized diversity
#3. effective species 
#====note
#====that the hierarchy for alpha diversity and effective species is almost completely identical

#==============================
#1. alpha diversity hierarchy 
#==== raw reads 
#bargraph
#wga vs non-wga
ggplot(t,aes(alpha.hierarchy,raw.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(alpha.hierarchy,raw.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== usable reads counts
ggplot(t,aes(alpha.hierarchy,usable.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(alpha.hierarchy,usable.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== unique seq
ggplot(t,aes(alpha.hierarchy,uniq_usable_seq.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(alpha.hierarchy,uniq_usable_seq.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== clonotype count
ggplot(t,aes(alpha.hierarchy,cln.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(alpha.hierarchy,cln.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#normelized diversity 
ggplot(t,aes(alpha.hierarchy,diversity_normed.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
ggplot(F.wga,aes(alpha.hierarchy,diversity_normed.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#effective species 
ggplot(t,aes(alpha.hierarchy,effective_species.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
ggplot(F.wga,aes(alpha.hierarchy,effective_species.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#===================================
#2. normalized diversity hierarchy
#==== raw reads 
#bargraph
#wga vs non-wga
ggplot(t,aes(normalized.hierarchy,raw.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(normalized.hierarchy,raw.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== usable reads counts
ggplot(t,aes(normalized.hierarchy,usable.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(normalized.hierarchy,usable.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== unique seq
ggplot(t,aes(normalized.hierarchy,uniq_usable_seq.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(normalized.hierarchy,uniq_usable_seq.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== clonotype count
ggplot(t,aes(normalized.hierarchy,cln.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(normalized.hierarchy,cln.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#effective species 
ggplot(t,aes(normalized.hierarchy,effective_species.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
ggplot(F.wga,aes(normalized.hierarchy,effective_species.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#===============================
#3. effective species  hierarchy
#==== raw reads 
#bargraph
#wga vs non-wga
ggplot(t,aes(eff_species.hierarchy,raw.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(eff_species.hierarchy,raw.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== usable reads counts
ggplot(t,aes(eff_species.hierarchy,usable.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(eff_species.hierarchy,usable.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== unique seq
ggplot(t,aes(eff_species.hierarchy,uniq_usable_seq.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(eff_species.hierarchy,uniq_usable_seq.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#==== clonotype count
ggplot(t,aes(eff_species.hierarchy,cln.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
#nonwga -> ca vs cf
ggplot(F.wga,aes(eff_species.hierarchy,cln.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))

#alpha diversity 
ggplot(t,aes(eff_species.hierarchy,alpha_diversity.count, fill=wga))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))
ggplot(F.wga,aes(eff_species.hierarchy,alpha_diversity.count, fill=sampleFraction))+geom_bar(position = "dodge",stat="identity")+theme(axis.text.x = element_text(angle = 90,hjust = 1))


#======================================================================================
#comparing all seq results (seq 21+24+26)
# no WGA 
#======================================================================================

# DNA yields
pdf(paste0(targetDir,'DNA yields (all seq).pdf'))
ggplot(F.wga,aes(sampleFraction,pre.dna,fill=seq.run))+facet_grid(.~submissionId)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#raw reads 
#boxplot with more then one layer of organization 
pdf(paste0(targetDir,'raw reads count (all seq).pdf'))
ggplot(F.wga,aes(pcrId,raw.count))+facet_grid(.~seq.run)+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#raw reads 

# correlation and summary 
summary(F.wga$raw.count)
summary(F.wga$usable.count)
summary(F.wga$cln.count)
summary(F.wga$alpha_diversity.count)
summary(F.wga$effective_species.count)
summary(F.wga$pred.lym.cyto)
summary(F.wga$hemo.cells.L)
summary(F.wga$pre.dna)

#ca vs cf 
#read
t.test (F.wga.ca$raw.count,F.wga.cf$raw.count)
t.test (F.wga.ca$usable.count,F.wga.cf$usable.count)
summary (F.wga.ca$raw.count)
summary (F.wga.cf$raw.count)
summary (F.wga.ca$usable.count)
summary (F.wga.cf$usable.count)
cor.test(F.wga.ca$raw.count,F.wga.cf$raw.count)
cor.test(F.wga.ca$usable.count,F.wga.cf$usable.count)
#clnotype
t.test (F.wga.ca$cln.count,F.wga.cf$cln.count)
summary(F.wga.ca$cln.count)
summary(F.wga.cf$cln.count)
#diversity
t.test (F.wga.ca$alpha_diversity.count,F.wga.cf$alpha_diversity.count)
summary(F.wga.ca$alpha_diversity.count)
summary(F.wga.cf$alpha_diversity.count)
t.test (F.wga.ca$effective_species.count,F.wga.cf$effective_species.count)
summary(F.wga.ca$effective_species.count)
summary(F.wga.cf$effective_species.count)

#wga -> make sure you take out run 1 form the data before running 
#raw reads
summary(T.wga$raw.count)
summary(F.wga$raw.count)
t.test(T.wga$raw.count,F.wga$raw.count)
#Usable reads 
t.test(T.wga$usable.count,F.wga$usable.count)
summary(F.wga$usable.count)
summary(T.wga$usable.count)
#clonotype
t.test(T.wga$cln.count,F.wga$cln.count)
summary(T.wga$cln.count)
summary(F.wga$cln.count)
#Alpha diversity
summary(T.wga$alpha_diversity.count)
summary(F.wga$alpha_diversity.count)
t.test(T.wga$alpha_diversity.count,F.wga$alpha_diversity.count)
#Effective species 
t.test (T.wga$effective_species.count,F.wga$effective_species.count)
summary(T.wga$effective_species.count)
summary(F.wga$effective_species.count)

#wga fraction analysis 
#raw reads
summary(T.wga.ca$raw.count)
summary(T.wga.cf$raw.count)
summary(F.wga.cf$raw.count)
summary(T.wga.cf$raw.count)
t.test (T.wga.ca$raw.count,F.wga.ca$raw.count)
t.test (T.wga.cf$raw.count,F.wga.cf$raw.count)
t.test (T.wga.ca$raw.count,T.wga.cf$raw.count)
#usable reads
summary(T.wga.ca$usable.count)
summary(T.wga.cf$usable.count)
summary(F.wga.ca$usable.count)
summary(F.wga.cf$usable.count)
t.test (T.wga.ca$usable.count,F.wga.ca$usable.count)
t.test (T.wga.cf$usable.count,F.wga.cf$usable.count)
t.test (T.wga.ca$usable.count,T.wga.cf$usable.count)
#clonotype 
summary(T.wga.ca$cln.count)
summary(T.wga.cf$cln.count)
summary(F.wga.cf$cln.count)
summary(F.wga.ca$cln.count)
t.test (T.wga.ca$cln.count,F.wga.ca$cln.count)
t.test (T.wga.cf$cln.count,F.wga.cf$cln.count)
t.test (T.wga.ca$cln.count,T.wga.cf$cln.count)
#diversity 
summary(F.wga.ca$alpha_diversity.count)
summary(F.wga.cf$alpha_diversity.count)
summary(T.wga.cf$alpha_diversity.count)
summary(T.wga.ca$alpha_diversity.count)
t.test (T.wga.ca$alpha_diversity.count,F.wga.ca$alpha_diversity.count)
t.test (T.wga.cf$alpha_diversity.count,F.wga.cf$alpha_diversity.count)
t.test (T.wga.ca$alpha_diversity.count,T.wga.cf$alpha_diversity.count)
#Effective species 
summary(F.wga.ca$effective_species.count)
summary(F.wga.cf$effective_species.count)
summary(T.wga.cf$effective_species.count)
summary(T.wga.ca$effective_species.count)
t.test (T.wga.ca$effective_species.count,F.wga.ca$effective_species.count)
t.test (T.wga.cf$effective_species.count,F.wga.cf$effective_species.count)
t.test (T.wga.ca$effective_species.count,T.wga.cf$effective_species.count)










