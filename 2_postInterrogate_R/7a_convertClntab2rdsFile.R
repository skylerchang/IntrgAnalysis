#this script converts clntab files into rds files for downstream processing in R
# -only important columns are kept -> file size reduction
# -depending on whether V and/or J genes could be identified, sequences are split into the following categories and printed to separate files:
# 1) 'vAndJ'
# 2) 'jOnly'
# 3) 'neitherVnorJ'
#
#the script also prints a summary statistics and plot for these categories


library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(seqinr)


targetFolder<-"../../Data/Clntab/"
outFolder<-"../../Results/"
outFolderDataRds<-"../../Data/Clntab_RDS/"
outFolderDataFasta<-"../../Data/Clntab_Fasta/"

datalist<-list()

#list all files in Clntab folder, which have long suffix
files<-list.files(targetFolder,pattern ="*.clntab",recursive = T)
#truncate suffix
files_short<-basename(sub("_L001_R1_001.fastq.processed.junctioned.profiled.clntab","",files))


plotlist<-list()

datalist_vAndJ<-list()
datalist_jOnly<-list()
datalist_neitherVnorJ<-list()

for (i in 1:length(files)){
#for (i in 2){
  print(files[i])
  t<-read_tsv(paste0(targetFolder,files[i]))
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

  d$junction<-factor(d$junction,levels=c('vAndJ','jOnly','neitherVnorJ'))
  
  #reduce number of lines to speed up plotting
  #dd<-d[sample(nrow(d),1000),]
  
  p<-ggplot(d,aes(junction,size,fill=vAndJchainSimplified))+geom_col()+scale_fill_brewer(palette='Set3')+theme(legend.position="top",legend.title=element_blank(),plot.title = element_text(hjust=0.5))+ggtitle(files_short[[i]])
  
  plotlist[[i]]<-ggplotGrob(p)
  
  datalist_vAndJ[[i]]<-d[d$junction=='vAndJ',]
  datalist_jOnly[[i]]<-d[d$junction=='jOnly',]
  datalist_neitherVnorJ[[i]]<-d[d$junction=='neitherVnorJ',]
}


#length analysis
for (i in 1:length(datalist_vAndJ)){
  t.vAndJ<-datalist_vAndJ[[i]]
  t.jOnly<-datalist_jOnly[[i]]
  t.neitherVnorJ<-datalist_neitherVnorJ[[i]]
  
  t.vAndJ$totalLength<-nchar(t.vAndJ$completeNtSeq)
  t.jOnly$totalLength<-nchar(t.jOnly$completeNtSeq)
  t.neitherVnorJ$totalLength<-nchar(t.neitherVnorJ$completeNtSeq)
  
  ggplot(t.vAndJ,aes(totalLength))+geom_bar()
  ggplot(t.jOnly,aes(totalLength))+geom_bar()
  ggplot(t.neitherVnorJ,aes(totalLength))+geom_bar()
}


#print rds files
saveRDS(list(datalist_vAndJ,files_short),paste0(outFolderDataRds,'clntab_vAndJ.rds'))
saveRDS(list(datalist_jOnly,files_short),paste0(outFolderDataRds,'clntab_jOnly.rds'))
saveRDS(list(datalist_neitherVnorJ,files_short),paste0(outFolderDataRds,'clntab_neitherVnorJ.rds'))

#print fasta files - vAndJ
dir.create(outFolderDataFasta)
dir.create(paste0(outFolderDataFasta,'/vAndJ'))
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

pdf(paste0(outFolder,'junctionVsNoju.pdf'))
marrangeGrob(plotlist,nrow=2,ncol=2)
dev.off()

seq<-c('aa','bb')
names<-c('name_aa','name_bb')
write.fasta(as.list(seq),names,paste0(outFolderDataRds,'/vAndJ.fasta'),as.string=F)

