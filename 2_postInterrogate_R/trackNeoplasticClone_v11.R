#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(vegan)
library(plyr)
library(here)
library(reshape2)

setwd(here())
getwd()

#=========== adjust =============

sampleNoCodesForFraction<-T        #T or F
thisIsReilly<-T

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ.rds'

loci<-c("IGH")
#===================================
t<-read_rds(rdsPath)

#datalist contains one tibble per sample
datalist<-t[[1]]
print(paste0("Number of files: ",length(datalist)))
#sample names are stored separately in a vector
files_short<-t[[2]]
files_short

indexSummary<-tibble()
templateSummary<-tibble()
diversity<-tibble()

for (i in 1:length(datalist)){
  a<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
  if(is.null(a)){next}
  if(nrow(a)==0){next}
  #============ standardize sample name =======
  #assume the following default format:
  # [1] "16-088089-3D1P3C1P1C1_Mcdermott-Gibbie"      
  # [2] "16-088089-3D1P4C1P1C1_Mcdermott-Gibbie"   
  #deviations:
  #[92] "k9-pc_D16-07-2FD4P59C1P1C1_OWN33-PAT1_S92"      
  #[93] "k9-pc-st_D16-07-2FD1P32C1P1C1_OWN33-PAT1_S93"  
  #<sampleId>_<sampleName>
  file<-sub("_S[0-9]+$","",files_short[i])
  #if this is reilly
  if (thisIsReilly==T){
    file<-sub("[0-9]+_","",file)
    file<-paste0(file,"_McCormick-Reilly")
  }
  #tranform samples to fit default format -> only one '_' which needs to be between sampleId & sampleName
  file<-sub("k9-pc-st_D16-07","D16-07-st",file)
  file<-sub("k9-pc_D16-07","D16-07",file)
  file<-sub("k9-nt-st_H2O","ntc",file)
  file<-sub("k9-nt_H2O","ntc-st",file)
  sampleId<-sub("_.*$","",file)
  sampleName<-sub("^.*_","",file)
  
  if (thisIsReilly==T){
    sampleId<-paste0(sampleId,"C1P1C1")
    file<-paste(sampleId,sampleName,sep="_")
  }
  #============ extract replicate & sample =======
  replicate<-ifelse(grepl("D[0-9]+P[0-9]*[02468]C",sampleId),"rep1","rep2")
  if(sampleNoCodesForFraction==T){
    fraction<-ifelse(grepl("-1D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$",sampleId),"PBMCs","Plasma")
  }else{
    fraction<-"fraction1"
  }
  sampleIdShort<-sub("D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$","",sampleId)
  submission<-sub("-[0-9a-zA-Z]+$","",sampleIdShort)
  
  #============ check for index sequences =======
  a$index<-'no'
  a$index[grepl("CARADYYDSFWAAFGYW",a$aaSeq)]<-'bella'
  a$index[grepl("CVTSVFSPW",a$aaSeq)]<-'bentley'
  a$index[grepl("CAPEIGWAAEYW",a$aaSeq)]<-'cheech'
  a$index[grepl("CAKDYYGTRNIYSLDYW",a$aaSeq)]<-'daisy'
  a$index[grepl("CVRSRRGYPETV#GMDYW",a$aaSeq)]<-'eddie1'
  a$index[grepl("CAKDGLVVPTGSIEYW",a$aaSeq)]<-'eddie2'
  a$index[grepl("CAKEWSSSNWYGDVGMDYW",a$aaSeq)]<-'gypsy'
  a$index[grepl("CARDYHGISWLEYW",a$aaSeq)]<-'jake'
  a$index[grepl("CVPFSPYGSWFADDQW",a$aaSeq)]<-'marishka'
  a$index[grepl("CVRGANSWFPDFDYW",a$aaSeq)]<-'reilly'
  a$index[grepl("CAKELYYDSYSVDYW",a$aaSeq)]<-'sona'
  
  #============ check for templates =======
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
  
  #split in separate tibbles
  a.withoutTemplate<-a[a$template=='no',]
  a.templateOnly<-a[a$template!='no',]
  
  #========== summarize index sequences ==============
  temp<-ddply(a.withoutTemplate[,c("size","index")],"index",na.rm=T, numcolwise(sum))
  if (nrow(temp)==0){next}
  temp$totalReads<-sum(temp$size)
  temp$percentage<-temp$size/temp$totalReads*100
  temp$submission<-submission
  temp$replicate<-replicate
  temp$fraction<-fraction
  indexSummary<-bind_rows(indexSummary,temp)
  
  #============ summarize templates for sample and store in tibble =======
  temp<-ddply(a.templateOnly,"template",numcolwise(sum))
  if(nrow(temp)==0){next}
  temp<-as_tibble(temp[,c("template","size")])
  temp$sample<-files_short[i]
  templateSummary<-bind_rows(templateSummary,temp)
  
  
  #============ diversity =======
  for (j in 1:length(loci)){
    aa<-as_tibble(ddply(a,"aaSeq",numcolwise(sum)))[,c("aaSeq","size")]
    aa.shannon<-diversity(aa$size, index = "shannon",base = exp(1))
    aa.simpson<-diversity(aa$size, index = "simpson")
    aa.invsimpson<-diversity(aa$size, index = "invsimpson")
    
    temp<-tibble(shannon=aa.shannon,simpson=aa.simpson,invsimpson=aa.invsimpson,sample=files_short[i],locus=loci[j],replicate=replicate,submission=submission,sampleIdShort=sampleIdShort)
    diversity<-bind_rows(diversity,temp)
  }
}

#======== diversity summary ===========
d<-diversity[rowSums(is.na(diversity)),]
#transform to long format
d.long<-as_tibble(melt(d,id.vars=c("sample","locus","replicate","submission","sampleIdShort")))

ggplot(d.long,aes(variable,value))+geom_boxplot()
d.long$value
#======== template summary ===========
templateSummary
ggplot(templateSummary,aes(template,size,fill=sample))+geom_col()+scale_fill_brewer(palette="Paired")+coord_flip()
ggplot(templateSummary,aes(sample,size,fill=template))+geom_col()+scale_fill_brewer(palette="Dark2")+coord_flip()

#======== index summary ===========
i<-indexSummary
i.bella<-i[i$index=="bella" & grepl("bella",i$submission),]
i.bentley<-i[i$index=="bentley" & grepl("bentle",i$submission),]
i.cheech<-i[i$index=="cheech" & grepl("cheech",i$submission),]
i.daisy<-i[i$index=="daisy" & grepl("daisy",i$submission),]
i.eddie1<-i[i$index=="eddie1" & grepl("eddie",i$submission),]
i.eddie2<-i[i$index=="eddie2" & grepl("eddie",i$submission),]
i.gypsy<-i[i$index=="gypsy" & grepl("gypsy",i$submission),]
i.jake<-i[i$index=="jake" & grepl("jake",i$submission),]
i.marishka<-i[i$index=="marishka" & grepl("marish",i$submission),]
i.reilly<-i[i$index=="reilly" & grepl("reilly",i$submission),]
i.sona<-i[i$index=="sona" & grepl("sona",i$submission),]

#total reads
ggplot(i.bentley,aes(submission,totalReads,fill=replicate))+geom_col()+facet_wrap(~fraction,ncol=1)
ggplot(i.cheech,aes(submission,totalReads,fill=replicate))+geom_col()+facet_wrap(~fraction,ncol=1)
ggplot(i.gypsy,aes(submission,totalReads,fill=replicate))+geom_col()+facet_wrap(~fraction,ncol=1)

#bella
labels<-c(17,16,13,10,7,4,0)
limits<-paste0("bella-0",seq(0,6,1))
pdf("../Results/Clones/bella.pdf")
ggplot(i.bella,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
 dev.off()

#bentley
labels<-c(4,2,0)
limits<-paste0("bentle-0",seq(0,2,1))
pdf("../Results/Clones/bentley.pdf")
ggplot(i.bentley,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#cheech
labels<-c(14,12,9,7,4,0)
limits<-paste0("cheech-0",seq(0,5,1))
pdf("../Results/Clones/cheech.pdf")
ggplot(i.cheech,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#daisy
labels<-c(38,35,34,32,30,27,23,18,14,12,4,0)
limits<-c(paste0("daisy-0",seq(0,9,1)),paste0("daisy-",seq(10,11,1)))
pdf("../Results/Clones/daisy.pdf")
ggplot(i.daisy,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#eddie1
labels<-c(18,17,15,12,9,6,3,0)
limits<-paste0("eddie-0",seq(0,7,1))
pdf("../Results/Clones/eddie1.pdf")
ggplot(i.eddie1,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#eddie2
labels<-c(18,17,15,12,9,6,3,0)
limits<-paste0("eddie-0",seq(0,7,1))
pdf("../Results/Clones/eddie2.pdf")
ggplot(i.eddie2,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#gypsy
labels<-c(31,29,26,24,21,17,13,9,7,2,0)
limits<-c(paste0("gypsy-0",seq(0,9,1)),"gypsy-10")
pdf("../Results/Clones/gypsy.pdf")
ggplot(i.gypsy,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#jake
#13 timpoints in timeline sheet, only 12 sequenced -> no relapse sample; last sample is 1w before relapse
labels<-c(49,47,43,41,38,34,30,26,24,12,6,1,0)
limits<-c(paste0("jake-0",seq(0,9,1)),paste0("jake-",seq(10,12,1)))
pdf("../Results/Clones/jake.pdf")
ggplot(i.jake,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")+
  annotate(geom="text",x=13,y=15,label="No sample",angle=90)

#marishka
labels<-c(28,26,23,21,18,14,8,6,0)
limits<-c(paste0("marish-0",seq(0,8,1)))
pdf("../Results/Clones/marishka.pdf")
ggplot(i.marishka,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#reilly - including post relapse
#13 samples sequenced (00-12); samples to & including replapse: 10 => 3 additional time points post relapse
labels<-c(32,30,27,25,22,16,12,8,6,0,-2,-5,-7)
limits<-c(paste0("reilly-0",seq(0,9,1)),paste0("reilly-",seq(10,12,1)))
pdf("../Results/Clones/reilly_inclPostRelapse.pdf")
ggplot(i.reilly,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#reilly
#13 samples sequenced (00-12); samples to & including replapse: 10 => 3 additional time points post relapse
labels<-c(32,30,27,25,22,16,12,8,6,0)
limits<-paste0("reilly-0",seq(0,9,1))
pdf("../Results/Clones/reilly.pdf")
ggplot(i.reilly,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()

#sona
labels<-c(8,7,6,5,2,0)
limits<-paste0("sona-0",seq(0,5,1))
pdf("../Results/Clones/sona.pdf")
ggplot(i.sona,aes(submission,percentage,fill=replicate))+
  geom_col(position="dodge")+
  facet_wrap(~fraction,ncol=1)+
  scale_fill_brewer(palette="Paired")+
  scale_x_discrete(limits=limits,labels=labels)+
  xlab("Time till relapse [weeks]")+ylab("Reads of neoplastic clone [%]")
dev.off()



i<-i[i$index=="reilly",]
i<-subset(i,select=(-index))

#pbmc
i.pbmc<-i[i$fraction=="fraction1",]
colnames(i.pbmc)<-paste(colnames(i.pbmc),"pbmc",sep=".")
i.pbmc<-subset(i.pbmc,select=(-fraction.pbmc))

i.pbmc.rep1<-i.pbmc[i.pbmc$replicate.pbmc=='rep1',]
colnames(i.pbmc.rep1)<-paste(colnames(i.pbmc.rep1),"rep1",sep=".")
i.pbmc.rep1<-subset(i.pbmc.rep1,select=(-replicate.pbmc.rep1))
colnames(i.pbmc.rep1)[4]<-"submission"
i.pbmc.rep1$submission<-as.factor(i.pbmc.rep1$submission)

i.pbmc.rep2<-i.pbmc[i.pbmc$replicate.pbmc=='rep2',]
colnames(i.pbmc.rep2)<-paste(colnames(i.pbmc.rep2),"rep2",sep=".")
i.pbmc.rep2<-subset(i.pbmc.rep2,select=(-replicate.pbmc.rep2))
colnames(i.pbmc.rep2)[4]<-"submission"
i.pbmc.rep2$submission<-as.factor(i.pbmc.rep2$submission)

#plasma
i.plasma<-i[i$fraction=="fraction2",]
colnames(i.plasma)<-paste(colnames(i.plasma),"plasma",sep=".")
i.plasma<-subset(i.plasma,select=(-fraction.plasma))

i.plasma.rep1<-i.plasma[i.plasma$replicate.plasma=='rep1',]
colnames(i.plasma.rep1)<-paste(colnames(i.plasma.rep1),"rep1",sep=".")
i.plasma.rep1<-subset(i.plasma.rep1,select=(-replicate.plasma.rep1))
colnames(i.plasma.rep1)[4]<-"submission"
i.plasma.rep1$submission<-as.factor(i.plasma.rep1$submission)

i.plasma.rep2<-i.plasma[i.plasma$replicate.plasma=='rep2',]
colnames(i.plasma.rep2)<-paste(colnames(i.plasma.rep2),"rep2",sep=".")
i.plasma.rep2<-subset(i.plasma.rep2,select=(-replicate.plasma.rep2))
colnames(i.plasma.rep2)[4]<-"submission"
i.plasma.rep2$submission<-as.factor(i.plasma.rep2$submission)

pbmc<-merge(i.pbmc.rep1,i.pbmc.rep2,by="submission")
plasma<-merge(i.plasma.rep1,i.plasma.rep2,by="submission")
i.rds<-as_tibble(merge(pbmc,plasma,by="submission"))

path<-"/Users/SKeller/Documents/Projects/MRD/DataAnalysis/MRD/Tables/ClonesInBlood/test.rds"
saveRDS(i.rds,path)

i.rds.long<-melt(i.rds,id.vars="submission")
i.rds.long$replicate<-ifelse(grepl("rep1",i.rds.long$variable),"rep1","rep2")
i.rds.long$fraction<-ifelse(grepl("pbmc",i.rds.long$variable),"pbmc","plasma")
i.rds.long$percentage<-ifelse(grepl("percentage",i.rds.long$variable),T,F)
ggplot(i.rds.long[i.rds.long$percentage==T,],aes(submission,value,fill=replicate))+geom_col()

       