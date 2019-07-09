#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(vegan)
library(plyr)
library(here)
library(reshape2)
library(openxlsx)
library(RColorBrewer)
library(ggforce)

setwd(here())
getwd()

#=========== adjust =============
run<-30
sampleNoCodesForFraction<-F        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'
resultsPath<-'../Results/Diversity/'

dir.create(resultsPath)

loci<-c("IGH","TRB","IGH+TRB")
#===================================
t<-read_rds(rdsPath)

#datalist contains one tibble per sample
datalist<-t[[1]]
length(datalist)
#sample names are stored separately in a vector
files_short<-t[[2]]
length(files_short)

#filter out WGA samples
includedLogical<-!grepl("WGA",files_short)
files_short<-files_short[includedLogical]
datalist<-datalist[includedLogical]

#initialize tibbles
templateSummary<-tibble()
diversity<-tibble()
readCount<-tibble()
#tibbles to aggregate the data from each loop
data.woTemp<-tibble()
data.tempOnly<-tibble()

for (i in 1:length(files_short)){  
  print(files_short[i])
  a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
  total<-nrow(a)
  if(nrow(a)==0){next}
  #============ standardize sample name =======
  #assume the following default format:
  # [1] "16-088089-3D1P3C1P1C1_Mcdermott-Gibbie"      
  # [2] "16-088089-3D1P4C1P1C1_Mcdermott-Gibbie"   
  #deviations:
  #[92] "k9-pc_D16-07-2FD4P59C1P1C1_OWN33-PAT1_S92"      
  #[93] "k9-pc-st_D16-07-2FD1P32C1P1C1_OWN33-PAT1_S93"  
  #<sampleId>_<sampleName>
  #tranform samples to fit default format -> only one '_' which needs to be between sampleId & sampleName
  file<-sub("_S[0-9]+$","",files_short[i])
  file<-sub("k9-pc-st_D16-07","D16-07-st",file)
  file<-sub("k9-pc_D16-07","D16-07",file)
  file<-sub("k9-nt-st_H2O","ntc",file)
  file<-sub("k9-nt_H2O","ntc-st",file)
  a$sampleId<-sub("_.*$","",file)
  a$sampleName<-sub("^.*_","",file)
  
  #============ extract replicate & sample =======
  a$replicate<-ifelse(grepl("D[0-9]+P[0-9]*[13579]",a$sampleId),"rep1","rep2")
  if(sampleNoCodesForFraction==T){
    a$fraction<-ifelse(grepl("-1D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$",a$sampleId),"fraction1","fraction2")
  }else{
    a$fraction<-"fraction1"
  }
  a$sampleIdShort<-sub("D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$","",a$sampleId)
  a$submission<-sub("-[0-9a-zA-Z]+$","",a$sampleIdShort)
  
  #============ filter seqs by seqs & ^C.{3,30}[FW]$ =======
  out<-nrow(a[!grepl("^C.{3,30}[FW]$",a$aaSeq),])
  perc<-out/total*100
  print(paste0("Filtered reads by '^C.{3,30}[FW]$': ",out,"/",total," (",perc,"%)"))
  a<-a[grepl("^C.{3,30}[FW]$",a$aaSeq),]
  if(nrow(a)==0){next}
  #============ split dataset into +/- templates =======
  a$template<-'no'
  #Akash‘s templates
  a$template[grepl("TGTGCATCACGACACAGTGGTCTGG",a$completeNtSeq)]<-'t1'
  a$template[grepl("TGTGCATCACGACCAGATCCACAGATCCATTGGTTACTGG",a$completeNtSeq)]<-'t2'
  #Tamara‘s templates
  a$template[grepl("TGTTCGCCTTATCGCCTTATGG",a$completeNtSeq)]<-'IGHV1-30_IGHJ4'
  a$template[grepl("TGTCTAGTACGCCTCTCTGCCTCTCTGCTAGTACGTGG",a$completeNtSeq)]<-'IGHV1-30_IGHJ6'
  a$template[grepl("TGTGCTTCTGCCTTTCTGCCTGCTCAGGATTCTGCCTGCTCAGGAGCTCAGGATTCTCTGG",a$completeNtSeq)]<-'IGHV3-1_IGHJ4'
  a$template[grepl("TGTGAGGAGTCCGTAGAGAGAGGAGTCCAGCGTAGCCATGCCTAAGGAGTCCCAGCCTCGGTAGAGAGAGCGCTGG",a$completeNtSeq)]<-'IGHV3-1_IGHJ6'
  a$template[grepl("TGTAAGGAGTAACTGCATAACTGCATACTAAGCCTAAGGAGTATGG",a$completeNtSeq)]<-'IGHV4-1_IGHJ4'
  a$template[grepl("TGTGTGCCTCTTTCCTCTACTAGATCGCCTCTCTATTATCCTCTAGAGTAGAGTAAGGAGTAGATCGCTATCCTCTGTAAGGAGTCCTCTACCTGG",a$completeNtSeq)]<-'IGHV4-1_IGHJ6'
  
  #split a into w/o template and template only
  a.woTemp<-a[a$template=='no',]
  a.tempOnly<-a[a$template!='no',]
  
  data.woTemp<-bind_rows(data.woTemp,a.woTemp)
  data.tempOnly<-bind_rows(data.tempOnly,a.tempOnly)
  
  template<-nrow(a.tempOnly)
  noTemplate<-nrow(a.woTemp)
  
  temp.count<-tibble(sample=a$sampleIdShort[1],total=total,CWFilterOut=out,template=template,noTemplate=noTemplate)
  readCount<-bind_rows(readCount,temp.count)
  
  #============ summarize templates for sample and store in tibble =======
  temp<-as_tibble(ddply(a.tempOnly[,c("template","size")],"template",numcolwise(sum)))
  if(nrow(temp)!=0){
    temp$sample<-files_short[i]
    templateSummary<-bind_rows(templateSummary,temp)
  }
  
  #============ diversity =======
  for (j in 1:length(loci)){
    #subset by locus if locus is NOT "IGH+TRB"
    if (loci[j]=="IGH+TRB"){
      aa<-a.woTemp
    }else{
      aa<-a.woTemp[a.woTemp$vAndJchainSimplified==loci[j],]
      }
    aa<-as_tibble(ddply(aa[,c("aaSeq","size")],"aaSeq",numcolwise(sum)))
    shannon<-diversity(aa$size,index = "shannon",base = exp(1))
    simpson<-diversity(aa$size,index = "simpson")
    invsimpson<-diversity(aa$size,index = "invsimpson")
    totalSizeByLocus<-sum(aa$size)
    clonotypeCountByLocus<-nrow(aa)
    effectiveSpecies<-exp(shannon)
    
    temp<-tibble(shannon=shannon,simpson=simpson,invsimpson=invsimpson,sample=files_short[i],locus=loci[j],replicate=a$replicate[1],submission=a$submission[1],sampleIdShort=a$sampleIdShort[1],size=totalSizeByLocus,clonotypeCount=clonotypeCountByLocus,effectiveSpecies=effectiveSpecies)
    diversity<-bind_rows(diversity,temp)
  }
}
#============ read count summary ===========
#print xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",readCount)
saveWorkbook(wb,"../Results/ReadCounts/readCountSummary.xlsx",overwrite=T)

readCount.long<-melt(readCount[,c(1,3:5)],id.vars="sample")
pdf('../Results/ReadCounts/readCountSummary.pdf')
ggplot(readCount.long,aes(sample,value,fill=variable))+geom_col()+coord_flip()
dev.off()

#============ locus stats ==========
wot<-data.woTemp[data.woTemp$vAndJchainSimplified=='IGH' | data.woTemp$vAndJchainSimplified=='TRB',]

wot.summary<-as_tibble(ddply(wot,c("vAndJchainSimplified","submission","replicate","fraction"),numcolwise(sum)))
wot.summary<-subset(wot.summary,select=-aaLength)

pdf("../Results/LocusStats/locusSummary_bySample.pdf")
ggplot(wot.summary,aes(vAndJchainSimplified,size,fill=replicate))+geom_col(position=position_dodge(preserve="single"))+facet_wrap(~submission)
dev.off()

if (run==28){
  wot.test<-wot.summary[-c(47,47,59,60,115,116,127,128),]
  pdf("../Results/LocusStats/locusSummary(19-018077,18-095326 removed).pdf")
  ggplot(wot.test,aes(vAndJchainSimplified,size,fill=replicate))+geom_col(position=position_dodge(preserve="single"))+facet_wrap(~submission)
  dev.off()
}

#convert to wide format - define function
transformToWide<-function(x){
  wot.summary.sub<-wot.summary[wot.summary$vAndJchainSimplified==x,]
  wot.summary.sub<-subset(wot.summary.sub,select=-vAndJchainSimplified)
  dcast(wot.summary.sub,submission+fraction~replicate,value.var="size")
  #colnames(wot.wide)[3:4]<-c(paste0(x,".rep1"),paste0(x,".rep2"))
}

#convert to wide format - use function
igh<-transformToWide("IGH")
trb<-transformToWide("TRB")
c<-merge(igh,trb,by=c("submission","fraction"))
colnames(c)[3:6]<-c("IGH.rep1","IGH.rep2","TRB.rep1","TRB.rep2")

#calculate replicate means
c$mean.IGH<-(c$IGH.rep1+c$IGH.rep2)/2
c$mean.TRB<-(c$TRB.rep1+c$TRB.rep2)/2

#calculate the igh/trb ratio 
c$locusRatio<-c$mean.IGH/c$mean.TRB

#plot ratios against each other 
#no log
pdf("../Results/LocusStats/locusRatio.pdf")
ggplot(c,aes(locusRatio,submission))+geom_point()
dev.off()
#log axis 
pdf("../Results/LocusStats/locusRatio-logAxis.pdf")
ggplot(c,aes(locusRatio,submission))+geom_point()+scale_x_log10()+ geom_vline(xintercept = 1)
dev.off()
#log transformation
pdf("../Results/LocusStats/locusRatio-logTransform.pdf")
ggplot(c,aes(log(locusRatio),submission))+geom_point()+ geom_vline(xintercept = 0)
dev.off()

#summary
summary(c$mean.IGH)
summary(c$mean.TRB)
summary(c$locus.Ratio)

#save as rds
saveRDS(c,"../Results/LocusStats/locusSummary.rds")

#print to xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",c)
saveWorkbook(wb,"../Results/LocusStats/locusSummary.xlsx",overwrite = T)

#============= diversity summary ===========
d<-diversity#[rowSums(is.na(diversity)),]
d.split<-split(d,d$submission)
str(d.split)

#print to xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",d)
addWorksheet(wb,"summary_IGH+TRB_only")
d.igh.trb<-d[d$locus=="IGH+TRB",]
writeData(wb,"summary_IGH+TRB_only",d.igh.trb)
for (i in 1:length(d.split)){
  addWorksheet(wb,d.split[[i]]$submission[1])
  writeData(wb,d.split[[i]]$submission[1],d.split[i],colNames = TRUE)
}
saveWorkbook(wb,paste0(resultsPath,"diversitySummary.xlsx"),overwrite = T)

#transform to long format
d.long<-as_tibble(melt(as.data.frame(d),id.vars=c("sample","locus","replicate","submission","sampleIdShort","size","clonotypeCount")))
d.long<-na.omit(d.long)
d.long.subset<-d.long[d.long$value!=Inf,]
d.long.subset<-d.long.subset[d.long.subset$locus!="IGH+TRB",]
d.long.subset$owner.patient<-sub("_S[0-9]+$","",d.long.subset$sample)
d.long.subset$owner.patient<-sub(".+_","",d.long.subset$owner.patient)
d.long.subset$owner<-sub("-.*","",d.long.subset$owner.patient)
str(d.long.subset)
d.long.subset$owner<-as.factor(d.long.subset$owner)

colorCount<-nlevels(d.long.subset$owner)
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))

#jitter all indexes
pdf(paste0(resultsPath,"diversityPlot_allIndexes.pdf"))
ggplot(d.long.subset,aes(variable,value))+
  geom_jitter(aes(shape=locus,color=owner),size=3)+
  facet_wrap(~variable, scales="free")+
  scale_color_manual(values=getPalette(colorCount))
dev.off()

#point - by size - all indexes 
pdf(paste0(resultsPath,"diversityPlot_allIndexes_withSize.pdf"))
ggplot(d.long.subset,aes(size,value))+
  geom_point(mapping=aes(shape=locus,color=owner),size=3)+
  facet_wrap(~variable,ncol=1,scales="free")+
  scale_color_manual(values=getPalette(colorCount))
dev.off()

#point - by size - shannon - with label
pdf(paste0(resultsPath,"diversityPlot_shannon_withSize_withLabel.pdf"))
d.long.subset2<-d.long.subset[d.long.subset$variable=="shannon",]
ggplot(d.long.subset2,aes(size,value))+
  geom_point(mapping=aes(shape=locus,color=owner),size=3)+
  facet_wrap(~variable,ncol=1,scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  geom_text(aes(label=owner),hjust=0, vjust=0,size=3)+
  scale_x_log10()+
  theme(legend.position="none")
dev.off

#point - by size - shannon - facet by owner
p<-ggplot(d.long.subset2,aes(size,value))+
  geom_point(mapping=aes(shape=locus,color=owner),size=3)+
  facet_wrap_paginate(~owner,ncol=3,nrow=4)+
  scale_color_manual(values=getPalette(colorCount))+
  scale_x_log10()+
  theme(legend.position="none")
pages<-n_pages(p)
pdf(paste0(resultsPath,"diversityPlot_shannon_withSize_facetByOwner.pdf"))
for (i in 1:pages){
  print(ggplot(d.long.subset2,aes(size,value))+
          geom_point(mapping=aes(shape=locus,color=owner),size=3)+
          facet_wrap_paginate(~owner,ncol=3,nrow=4,page=i)+
          scale_color_manual(values=getPalette(colorCount))+
          scale_x_log10()+
          theme(legend.position="none"))
}
dev.off()

#point - by size - simpson - with label
d.long.subset2<-d.long.subset[d.long.subset$variable=="simpson",]
pdf(paste0(resultsPath,"diversityPlot_simpson_withSize_withLabel.pdf"))
ggplot(d.long.subset2,aes(size,value))+
  geom_point(mapping=aes(shape=locus,color=owner),size=3)+
  facet_wrap(~variable,ncol=1,scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  geom_text(aes(label=owner),hjust=0, vjust=0,size=3)+
  scale_x_log10()+
  theme(legend.position="none")
dev.off()

#point - by size - invsimpson - with label
d.long.subset2<-d.long.subset[d.long.subset$variable=="invsimpson",]
pdf(paste0(resultsPath,"diversityPlot_invsimpson_withSize_withLabel.pdf"))
ggplot(d.long.subset2,aes(size,value))+
  geom_point(mapping=aes(shape=locus,color=owner),size=3)+
  facet_wrap(~variable,ncol=1,scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  geom_text(aes(label=owner),hjust=0, vjust=0,size=3)+
  scale_x_log10()+
  theme(legend.position="none")
dev.off()

#point - by size - effectiveSpecies - with label
d.long.subset2<-d.long.subset[d.long.subset$variable=="effectiveSpecies",]
pdf(paste0(resultsPath,"diversityPlot_effectiveSpecies_withSize_withLabel.pdf"))
ggplot(d.long.subset2,aes(size,value))+
  geom_point(mapping=aes(shape=locus,color=owner),size=3)+
  facet_wrap(~variable,ncol=1,scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  geom_text(aes(label=owner),hjust=0, vjust=0,size=3)+
  scale_x_log10()+
  theme(legend.position="none")
dev.off()


#======== template summary ===========
templateSummary
pdf('../Results/Templates/templatesBytemplate.pdf')
ggplot(templateSummary,aes(template,size,fill=sample))+geom_col()+coord_flip()
dev.off()
pdf('../Results/Templates/templatesBySample.pdf')
ggplot(templateSummary,aes(sample,size,fill=template))+geom_col()+coord_flip()
dev.off()
