library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)

targetFolder<-"../../Data/Clntab/"
outFolder<-"../../Data/CountPlots/"

datalist<-list()

files<-list.files(targetFolder,pattern ="*.clntab",recursive = T)
files_short<-basename(sub("_L001_R1_001.fastq.processed.junctioned.profiled.clntab","",files))

seqs<-tibble(aaSeq=character(),file=character())

for (i in 1:length(files)){
  print(files[i])
  t<-read_tsv(paste0(targetFolder,files[i]))
  tt<-subset(t, select=c("sequence.5-GENE","sequence.3-GENE","sequence.JUNCTION.aa seq","sequence.JUNCTION.aa seq.len","sequence.chain"))
  colnames(tt)<-c("v","j","aaSeq","aaLength","chain")
  conditionOne<-grepl("^C.{4,34}[WF]$",tt$aaSeq)
  conditionTwo<-!grepl("[\\*#]",tt$aaSeq)
  tt$validCdr<-ifelse(conditionOne & conditionTwo,TRUE,FALSE)
  tt$chain<-factor(tt$chain,levels=c('A','B','D','G','H','K','NA')) 
  
  #bind valid cdrs to 'seqs'
  temp2<-tt[tt$validCdr==TRUE,"aaSeq"]
  temp1<-tibble(file=rep(files_short[i],nrow(temp2)))
  temp<-cbind(temp1,temp2)

  seqs<-rbind(seqs,temp)
  datalist[[i]]<-tt
}

seqs<-as_tibble(seqs)
require(plyr)
seqs$concat<-paste(seqs$file,seqs$aaSeq)
summary<-ddply(seqs,.(aaSeq),summarise,freq=length(aaSeq))
summary<-head(summary[order(-summary$freq),],n=10)
topTen<-summary$aaSeq
seqs$aaSeq<-factor(seqs$aaSeq)
seqs$file<-factor(seqs$file)

seqsTop<-seqs[seqs$aaSeq %in% topTen,]


colorCount<-nlevels(seqs$file)
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))
pdf(paste0(outFolder,"topTenClones_cloneVsAbundance.pdf"))
ggplot(seqsTop,aes(aaSeq,fill=file))+geom_bar()+coord_flip()+scale_fill_manual(values = getPalette(colorCount))+guides(fill=guide_legend(title = 'Samples',ncol = 1))
dev.off()

pdf(paste0(outFolder,"topTenClones_sampleVsClone.pdf"))
ggplot(seqsTop,aes(file,aaSeq))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



plotlistLoci<-list()
plotlistTRB<-list()
plotlistIGH<-list()

for (i in 1:length(datalist)){
  t<-datalist[[i]]
  p<-ggplot(t,aes(chain, fill=(validCdr)))+geom_bar()+ggtitle(files_short[i])+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('A','B','D','G','H','K','NA'))
  plotlistLoci[[i]]<-ggplotGrob(p)
  
  t.b<-t[t$validCdr==T & t$chain=='B',]
  p<-plotlistTRB[[i]]<-ggplot(t.b,aes(aaLength,fill=(aaSeq)))+geom_bar()+guides(fill=F)+xlim(4,38)+ggtitle(files_short[i])+theme(plot.title = element_text(hjust = 0.5,size=10))
  plotlistTRB[[i]]<-ggplotGrob(p)
  
  t.h<-t[t$validCdr==T & t$chain=='H',]
  p<-plotlistIGH[[i]]<-ggplot(t.h,aes(aaLength,fill=(aaSeq)))+geom_bar()+guides(fill=F)+xlim(4,38)+ggtitle(files_short[i])+theme(plot.title = element_text(hjust = 0.5,size=10))
  plotlistIGH[[i]]<-ggplotGrob(p)
}


require(gridExtra)
pdf(paste0(outFolder,'clntab_byLocus.pdf'))
marrangeGrob(plotlistLoci, nrow=2, ncol=2)
dev.off()

pdf(paste0(outFolder,'clntab_trbByAaLength.pdf'))
print(marrangeGrob(plotlistTRB, nrow=4, ncol=4))
dev.off()

pdf(paste0(outFolder,'clntab_ighByAaLength.pdf'))
print(marrangeGrob(plotlistIGH, nrow=4, ncol=4))
dev.off()

ggplot(seqs,aes(aaSeq))+geom_bar()


#extract nonsense reads
i<-1
t<-read_tsv(paste0(targetFolder,files[i]))
tt<-subset(t, select=c("sequence.5-GENE","sequence.3-GENE","sequence.JUNCTION.aa seq","sequence.JUNCTION.aa seq.len","sequence.chain","sequence.nt seq"))
str(tt)  
colnames(tt)<-c("v","j","aaSeq","aaLength","chain","nt")
trbj16<-tt$nt[tt$j=='TRBJ1-6']
noju<-tt$nt[tt$aaSeq=='noju']
write.fasta(as.list(noju[1:100]),rep("trbj16",100),'trbj16_100.fa')
