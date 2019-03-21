library(tidyverse)
library(gridExtra)


t<-read_rds('../../Data/clntab.rds')
datalist<-t[[1]]
files_short<-t[[2]]

outpath<-'../../Results/Histograms/'

plotlist_histogramFacet<-list()
plotlist_density<-list()

for (i in 1:length(datalist)){
  s<-datalist[[i]]
  
  loci<-c("TRA","TRB","TRD","TRG")
  
  #eliminate reads with too long aaSeq
  s<-s[s$aaLength<=40,]
  #eliminate NA reads
  s<-s[rowSums(is.na(s))==0,]

  p<-ggplot(s,aes(aaLength,fill=vAndJchainSimplified))+geom_bar()+facet_wrap(~vAndJchainSimplified)+ggtitle(files_short[i])+theme(plot.title=element_text(hjust=0.5))+guides(fill=F)
  plotlist_histogramFacet[[i]]<-ggplotGrob(p)
  
#  p<-ggplot(s,aes(aaLength,fill=vAndJchainSimplified))+geom_density(alpha = 0.5, position = "stack",adjust=3)
#  plotlist_density[[i]]<-ggplotGrob(p)

}

pdf(paste0(outpath,'cdr3Length_histogram_all.pdf'))
marrangeGrob(plotlist_histogramFacet,nrow=1,ncol=1)
dev.off()

pdf(paste0(outpath,'cdr3Length_histogram_allOnOnePage.pdf'))
marrangeGrob(plotlist_histogramFacet,nrow=4,ncol=3)
dev.off()

pdf(paste0(outpath,'cdr3Length_density_all.pdf'))
marrangeGrob(plotlist_density,nrow=4,ncol=3)
dev.off()

#below is garbage

t<-datalist[[1]]

t.trg<-t[t$vAndJchainSimplified=='TRG',]

ggplot(t.trg,aes(aaLength,size))+geom_col()+facet_wrap(vGene~jGene)



t <- read.table("1A_S23/1A_S23_TRB_col15-21-23-27.clntab",header=F,fill = TRUE)
colnames(t)<-c("aaSeq","length","functionality","size")

t$aaSeq <- reorder(t$aaSeq, t$size)
t<-head(t, n=100)

reds<-c(hsv(.01,.9,.9),hsv(.01,.9,.85),hsv(.01,.9,.8),hsv(.01,.9,.75),hsv(.01,.9,.7))
greens<-c(hsv(.4,.9,.9),hsv(.4,.9,.85),hsv(.4,.9,.8),hsv(.4,.9,.75),hsv(.4,.9,.7))
t$color <- ifelse(t$functionality == "productive", sample(greens), sample(reds))

ggplot(data = t, aes(x = length, y = size, fill = aaSeq)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=t$color, guide = guide_legend(reverse=TRUE, limits='CSVLLCQLPGGR#YERYF')) 



ggplot(data = t, aes(x = length, y = size, fill = aaSeq)) + 
  geom_bar(stat = "identity") + scale_fill_brewer(palette="Dark2", guide = guide_legend(reverse=TRUE))

ggplot(data = t, aes(x = length, y = size, fill = aaSeq)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=t$color, guide = guide_legend(reverse=TRUE)) 











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
