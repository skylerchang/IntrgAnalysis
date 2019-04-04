#this script is specific to Run02; save under new name and adapt to other runs)
#I threw various pieces together; the code still needs to be adjusted

library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)

#========== adjust the following variables ================
t<-read_rds('../../Data/clntab_RDS/clntab_vAndJ.rds')
outpath<-'../../Results/Plots/'
loci<-c("TRA","TRB","TRD","TRG")

#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

#plotlist_histogramFacet<-list()
#plotlist_density<-list()

library(plyr)



for (i in 1:length(datalist)){
  t<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]

  nrow(t)
  #================ aaLength vs. aaSeq =================================
  #collapse reads with identical aaSeq
  library(plyr)
  for (j in 1:length(loci)){
    #subset for locus
    tt<-t[t$vAndJchainSimplified==loci[j],c("aaSeq","size")]
    tt<-tt[!is.na(tt$aaSeq),]
    #collapse lines with identical 'aaSeq' and sum up 'size'
    tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
    #sort in descending order of 'size'
    tt<-tt[order(-tt$size),]
    #convert the aaSeqs for all clones less abundant than 'n' into generic sequences
    n<-80
    tt$genericAaSeq<-as.vector(mapply(paste,mapply(rep,'x',nchar(tt$aaSeq)),collapse=''))
    #re-assign original aaSeq for first n clonotypes
    tt$aaSeq[(n+1):nrow(tt)]<-tt$genericAaSeq[(n+1):nrow(tt)]
    ttt<-ddply(tt,"aaSeq",numcolwise(sum))
    #assign color
    ttt$color<-"grey"
    b<-rep(brewer.pal(8,"Accent"),10)
    ttt$color[1:n]<-b[1:n]
    ttt$length<-nchar(ttt$aaSeq)
    
    #plot
    pdf(paste0(outpath,"/cdrAaSeqByAaLength_",names[i],"_size25_unicolor.pdf"))
    ggplot(ttt,aes(length,size,fill=aaSeq))+geom_col(position = position_stack(reverse = TRUE))+guides(fill=F)+scale_fill_manual(values=ttt$color)+labs(title=loci[j])+theme_bw(base_size = 25)
    dev.off()
  }
  
  #================ vGene usage =================================
  t
  tt<-as_tibble(ddply(t[,c("vGene","size")],"vGene",numcolwise(sum)))         
  tt<-tt[order(-tt$size),]
  tt$v<-factor(tt$v,levels=tt$v)
  
  #classify by subgroup
  tt$vSubgroup<-NA
  tt$vSubgroup[grepl("IGLV1",tt$v)]<-'V1'
  tt$vSubgroup[grepl("IGLV2",tt$v)]<-'V2'
  tt$vSubgroup[grepl("IGLV3",tt$v)]<-'V3'
  tt$vSubgroup[grepl("IGLV8",tt$v)]<-'V8'
  
  #plot all V genes in descending order of frequency
  pdf(paste0(outpath,"/vGeneUsage-individually_",names[i],".pdf"))
  ggplot(tt,aes(v,size,fill=vSubgroup))+geom_col()+scale_fill_brewer(palette = "Dark2")+labs(title=names[i])+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  #plot frequency of v gene subgroups
  pdf(paste0(outpath,"/vGeneUsage-bySubgroup_",names[i],".pdf"))
  colorCount<-nrow(tt)
  getPalette<-colorRampPalette(brewer.pal(9,'Set1'))
  ggplot(tt,aes(vSubgroup,size,fill=v))+geom_col()+guides(fill=F)+scale_fill_manual(values=sample(getPalette(colorCount)))+labs(title=names[i])+theme_bw(base_size = 15)
  dev.off()
  
  #================ jGene usage =================================
  t
  tt<-as_tibble(ddply(t[,c("j","size")],"j",numcolwise(sum)))         
  tt<-tt[order(-tt$size),]
  tt$j<-factor(tt$j,levels=tt$j)
  
  #plot all j genes in descending order of frequency
  pdf(paste0(outpath,"/jGeneUsage_",names[i],".pdf"))
  ggplot(tt,aes(j,size))+geom_col()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_brewer(palette = "Dark2")+labs(title=names[i],x="J gene usage", y="Sequence reads")+theme_bw(base_size = 15)
  dev.off()
}


#================== determine vGene rank ============================
data<-tibble()

#combine data from individual tibbles into one 'data' tibble
for (i in 1:length(datalist)){
  #subset tibble
  s<-datalist[[i]][,c("vGene",'size','vAndJchainSimplified')]
  #remove lines with ambiguous V genes or NA lines
  s<-s[!grepl('=',s$vGene),]
  s<-s[!is.na(s$vGene),]
  
  for (j in 1:length(loci)){

    t<-s[s$vAndJchainSimplified==loci[j],]
    
    #combine rows with same vGene
    tt<-ddply(t,"vGene",numcolwise(sum))
 
    tt<-tt[order(-tt$size),]
    tt$rank<-c(1:nrow(tt))
    tt$locus<-rep(loci[j],nrow(tt))
    tt$sample<-rep(files_short[i],nrow(tt))
    tt$organ<-NA
    tt$organ[grepl("L_",tt$sample)]<-'Lymph node'
    tt$organ[grepl("S_",tt$sample)]<-'Spleen'
    tt$organ[grepl("T_",tt$sample)]<-'Thymus'
    tt$individual<-sub("._.*","",tt$sample)
    
    #print to xlsx
    #writeData(wb, sheet = j,data)
    
    #limit to maximal top 10
    if (nrow(tt)<10){
      tt<-tt[1:nrow(tt),]
    }else{
      tt<-tt[1:10,]
    }
    data<-bind_rows(data,tt)
    }
}
detach("package:plyr",unload=T)

#print xlsx file to have the numerical data
wb <- createWorkbook()
addWorksheet(wb, 'test', gridLines = F)
writeData(wb, sheet = 1,data)
saveWorkbook(wb,paste0(outpath,'vGeneUsage_all.xlsx'),overwrite = TRUE)

#print V genes by rank
for (j in 1:length(loci)){
  d<-data[data$locus==loci[j],]
  d$organ<-as.factor(d$organ)
  
  #generate additinoal colors
  maxColors<-c(brewer.pal(8,"Accent"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"))
  colors<-maxColors[1:length(unique(d$vGene))]
  
  pdf(paste0(outpath,'vGenesByRank_',loci[j],'.pdf'))
  p<-ggplot(d,aes(individual,rank,group=vGene))+geom_line(aes(color=vGene))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_point(aes(color=vGene))+scale_color_manual(values = colors)+facet_wrap(~organ)+expand_limits(y=10)+scale_y_continuous(breaks=seq(10,1),trans = "reverse")+ggtitle(paste0("V gene usage by rank - ",loci[j]))+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

#===================== compare cdr3 length for different loci ===================
for (i in 1:length(datalist)){
  #subset tibble
  s<-datalist[[i]][,c("aaSeq",'size','vAndJchainSimplified')]
  #combine rows with same vGene
  tt<-ddply(t,"aaSeq",numcolwise(sum))
}
#adjust below piece
for (j in 1:length(loci)){
  s<-s[!is.na(s$vGene),]
  t<-s[s$vAndJchainSimplified==loci[j],]
  
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
