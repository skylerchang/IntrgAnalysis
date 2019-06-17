#this script converts clntab files into rds (& fasta files) for downstream processing in R
# -only important columns are kept -> file size reduction
# -depending on whether V and/or J genes could be identified, sequences are split into the following categories and printed to separate files:
# 1) 'vAndJ'
# 2) 'vOnly'
# 3) 'jOnly'
# 4) 'neitherVnorJ'
#
#the script also prints a summary statistics and plot for these categories


library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(seqinr)
library(ggforce)


targetFolder<-"../Data/Clntab/"
outFolder<-"../Results/ClntabJunctionCompare/"
outFolderDataRds<-"../Data/Clntab_RDS/"
outFolderDataFasta<-"../Data/Clntab_Fasta/"

#Dummy clntab file created for '18-069890-3D1P4C1P1C1' (run 28) since this there is no clntab file
#=>intestigate reason

setwd(here())
getwd()

dir.create(outFolder,recursive=T)

#list all files in Clntab folder, which have long suffix
files<-list.files(targetFolder,pattern ="*.clntab",recursive = T)
#truncate suffix
files_short<-basename(sub("_L001_R1_001.fastq.processed.junctioned.profiled.clntab","",files))

#============== validate file number & format of file names ==================
#the following assumes that samples are in duplicate and have the following format:
#[1] "16-088089-3D1P3C1P1C1_Mcdermott-Gibbie_S1"      
#[2] "16-088089-3D1P4C1P1C1_Mcdermott-Gibbie_S2"  
#check if number of files is even
if ((length(files_short) %% 2)!=0){stop("Number of files should be even (replicates expected). Exiting ...")}
#expand to identify replicates
#========================================================

#initialize things
datalist_vAndJ<-list()
datalist_vOnly<-list()
datalist_jOnly<-list()
datalist_neitherVnorJ<-list()
d.aggregated<-tibble()

for (i in 1:length(files)){
  print(files[i])
  t<-read_tsv(paste0(targetFolder,files[i]))
  
  #for files w/o reads -> insert empty tibble into datalists
  if(nrow(t)==0){
    df <- as_tibble(matrix(ncol = 7, nrow = 0))
    x <- c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")
    colnames(df) <- x
    datalist_vAndJ[[i]]<-df
    datalist_vOnly[[i]]<-df
    datalist_jOnly[[i]]<-df
    datalist_neitherVnorJ[[i]]<-df
    next
  }
  
  d<-subset(t, select=c("sequence.5-GENE","sequence.3-GENE","sequence.JUNCTION.aa seq","sequence.JUNCTION.aa seq.len","sequence.size","sequence.nt seq"))
  colnames(d)<-c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq")
  
  d$vChain<-NA
  d$vChain[grepl("IGH",d$vGene)]<-'IGH'
  d$vChain[grepl("TRA",d$vGene)]<-'TRA'
  d$vChain[grepl("TRB",d$vGene)]<-'TRB'
  d$vChain[grepl("TRD",d$vGene)]<-'TRD'
  d$vChain[grepl("TRG",d$vGene)]<-'TRG'
  
  d$jChain<-NA
  d$jChain[grepl("IGH",d$jGene)]<-'IGH'
  d$jChain[grepl("TRA",d$jGene)]<-'TRA'
  d$jChain[grepl("TRB",d$jGene)]<-'TRB'
  d$jChain[grepl("TRD",d$jGene)]<-'TRD'
  d$jChain[grepl("TRG",d$jGene)]<-'TRG'
  
  d$junction<-NA
  d$junction[d$vGene=='no5' & d$jGene=='no3']<-'neitherVnorJ'
  d$junction[d$vGene=='no5' & d$jGene!='no3']<-'jOnly'
  d$junction[d$vGene!='no5' & d$jGene=='no3']<-'vOnly'
  d$junction[d$vGene!='no5' & d$jGene!='no3']<-'vAndJ'
  
  d$vAndJchain<-paste(d$vChain,d$jChain,sep='-')
  
  d$vAndJchainSimplified<-NA
  d$vAndJchainSimplified[grepl('IGH',d$vAndJchain)]<-'IGH'
  d$vAndJchainSimplified[d$vAndJchain=='NA-NA' | d$vAndJchain=='TRD-TRA']<-'Other'
  d$vAndJchainSimplified[grepl('TRA',d$vAndJchain)]<-'TRA'
  d$vAndJchainSimplified[grepl('TRB',d$vAndJchain)]<-'TRB'
  d$vAndJchainSimplified[grepl('TRD',d$vAndJchain)]<-'TRD'
  d$vAndJchainSimplified[grepl('TRG',d$vAndJchain)]<-'TRG'

  d$junction<-factor(d$junction,levels=c('vAndJ','vOnly','jOnly','neitherVnorJ'))
  
  datalist_vAndJ[[i]]<-d[d$junction=='vAndJ',]
  datalist_vOnly[[i]]<-d[d$junction=='vOnly',]
  datalist_jOnly[[i]]<-d[d$junction=='jOnly',]
  datalist_neitherVnorJ[[i]]<-d[d$junction=='neitherVnorJ',]
  
  temp<-as_tibble(aggregate(d$size,by=list(junction=d$junction,locus=d$vAndJchainSimplified),FUN=sum))
  temp$sample<-files_short[i]
  d.aggregated<-bind_rows(d.aggregated,temp)
}

#================ plot junction types ===================
#shorten sample name; currently: "16-088089-3D1P3C1P1C1_Mcdermott-Gibbie_S1" 
unique(d.aggregated$sample)
d.aggregated$sample<-sub("_S[0-9]+$","",d.aggregated$sample)
d.aggregated$sample<-sub("_","\n",d.aggregated$sample)
d.aggregated$sample<-as.factor(d.aggregated$sample)

#all combined
pdf(paste0(outFolder,'junctionTypesCompared_all.pdf'))
ggplot(d.aggregated,aes(junction,x,fill=locus))+geom_col()
dev.off()

#separate as facet on one page
pdf(paste0(outFolder,'junctionTypesCompared_all_facet.pdf'))
ggplot(d.aggregated,aes(junction,x,fill=locus))+geom_col()+facet_wrap(d.aggregated$sample)
dev.off()

#separate as facet on multiple pages
p<-ggplot(d.aggregated,aes(junction,x,fill=locus))+geom_col()+facet_wrap_paginate(d.aggregated$sample,ncol=3,nrow=(4))
pages<-n_pages(p)
pdf(paste0(outFolder,'junctionTypesCompared_all_facet_3x4.pdf'))
for (i in 1:pages){
  print(ggplot(d.aggregated,aes(junction,x,fill=locus))+geom_col()+facet_wrap_paginate(d.aggregated$sample,ncol=3,nrow=(4),page=i))
}
dev.off()

#separate as facet on multiple pages; scales = free
p<-ggplot(d.aggregated,aes(junction,x,fill=locus))+geom_col()+facet_wrap_paginate(d.aggregated$sample,ncol=3,nrow=4)
pages<-n_pages(p)
pdf(paste0(outFolder,'junctionTypesCompared_all_facet_3x4_freeScales.pdf'))
for (i in 1:pages){
  print(ggplot(d.aggregated,aes(junction,x,fill=locus))+geom_col()+facet_wrap_paginate(d.aggregated$sample,ncol=3,nrow=4,page=i,scales="free"))
}
dev.off()

#point plot
d.aggregated$replicate<-ifelse(grepl("D[0-9]+P[0-9]*[02468]C",d.aggregated$sample),"rep1","rep2")
d.aggregated$submission<-sub("P[0-9]+C.*$","",d.aggregated$sample)
d.aggregated$owner.patient<-sub("_S[0-9]+$","",d.aggregated$sample)
d.aggregated$owner.patient<-sub("^.*_","",d.aggregated$owner.patient)
d.aggregated$id<-paste(d.aggregated$submission,d.aggregated$owner.patient,sep="\t")

p<-ggplot(d.aggregated,aes(replicate,x))+
  geom_point(aes(fill=locus,shape=junction))+
  scale_fill_brewer(palette="Set1")+
  coord_flip()+
  facet_wrap_paginate(~id,ncol=2,nrow=10,page=i)+
  scale_shape_manual(values=c(21,22,23,24))

pages<-n_pages(p)
pdf(paste0(outFolder,'junctionTypesCompared_all_facet_pointplot.pdf'))
for (i in 1:pages){
  print(ggplot(d.aggregated,aes(replicate,x))+
    geom_point(aes(fill=locus,shape=junction))+
    scale_fill_brewer(palette="Set1")+
    coord_flip()+
    facet_wrap_paginate(~id,ncol=2,nrow=10,page=i)+
    scale_shape_manual(values=c(21,22,23,24)))
}
dev.off()


#================ print rds files =====================
dir.create(outFolderDataRds)
saveRDS(list(datalist_vAndJ,files_short),paste0(outFolderDataRds,'clntab_vAndJ.rds'))
saveRDS(list(datalist_vOnly,files_short),paste0(outFolderDataRds,'clntab_vOnly.rds'))
saveRDS(list(datalist_jOnly,files_short),paste0(outFolderDataRds,'clntab_jOnly.rds'))
saveRDS(list(datalist_neitherVnorJ,files_short),paste0(outFolderDataRds,'clntab_neitherVnorJ.rds'))

#================ print fasta files =====================
#print fasta files - vAndJ
dir.create(paste0(outFolderDataFasta,'/vAndJ'),recursive=T)
for (i in 1:length(datalist_vAndJ)){
  header<-paste(datalist_vAndJ[[i]]$vGene,datalist_vAndJ[[i]]$jGene,sep='_')
  outFile<-paste0(outFolderDataFasta,'/vAndJ/',files_short[i],'_vAndJ.fasta')
  write.fasta(as.list(datalist_vAndJ[[i]]$completeNtSeq),header,outFile)
}

#print fasta files - neitherVnorJ
dir.create(paste0(outFolderDataFasta,'/neitherVnorJ'))
for (i in 1:length(datalist_neitherVnorJ)){
  header<-paste(datalist_neitherVnorJ[[i]]$vGene,datalist_neitherVnorJ[[i]]$jGene,sep='_')
  outFile<-paste0(outFolderDataFasta,'/neitherVnorJ/',files_short[i],'_neitherVnorJ.fasta')
  write.fasta(as.list(datalist_neitherVnorJ[[i]]$completeNtSeq),header,outFile)
}

