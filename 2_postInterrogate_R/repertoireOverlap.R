library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)
library(tcR)

#========== adjust the following variables ================
t<-read_rds('../../Data/clntab_RDS/clntab_vAndJ.rds')
outpath<-'../../Results/cdrAaLengthVsAaSeq'
loci<-c("TRA","TRB","TRD","TRG")
#Top n clones that are being displayed in separate colors (all other clones are grey)
n<-80

#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]
pcrId<-sub("_.*","",files_short)

data<-tibble()

for (i in 1:length(datalist)){
  t<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]
  t$sample<-pcrId[i]
  data<-bind_rows(data,t)
}

data<-data[data$aaLength<30,]
data<-data[grepl("^C.*[FW]$",data$aaSeq),]
data<-data[!data$vAndJchainSimplified=='TRA',]
data.igh<-data[data$vAndJchainSimplified=='IGH',]
data.trb<-data[data$vAndJchainSimplified=='TRB',]
ggplot(data,aes(aaLength,size))+geom_col()+facet_wrap(~vAndJchainSimplified)

#============ igh =================
pdf(paste0(outpath,'facet-igh.pdf'))
ggplot(data.igh,aes(aaLength,size))+geom_col()+facet_wrap(~sample)
dev.off()

s<-data[data$aaSeq=='CGCSWSATLEYW',]
ggplot(s,aes(sample,size))+geom_col()+coord_flip()
library(dplyr)
ddply(s,"sample",numcolwise(sum))

#============ trb =================
#cdr3 length for all samples
pdf(paste0(outpath,'facet-trb.pdf'))
ggplot(data.trb,aes(aaLength,size))+geom_col()+facet_wrap(~sample)
dev.off()

#sample '15-094420-2D' seems clonal and has a disproportionate number of reads <- get clone
data.trb.094420<-data.trb[grepl('15-094420-2D',data.trb$sample),]
data.trb.094420<-ddply(data.trb.094420,"aaSeq",numcolwise(sum))
head(data.trb.094420[order(-data.trb.094420$size),c('aaSeq','size')])
#check if neoplastic clone 'CASRVGGRNTLHF' appears in other samples
CASRVGGRNTLHF<-data.trb[data.trb$aaSeq=='CASRVGGRNTLHF',]
ddply(CASRVGGRNTLHF,"sample",numcolwise(sum))

#remove 15-094420-2D from pool 
data.trb.subset<-data.trb[!grepl('15-094420-2D',data.trb$sample),]
#plot juLength vs. juAaSeq
pdf(paste0(outpath,'facet-trbMinus15-094420.pdf'))
ggplot(data.trb.subset,aes(aaLength,size))+geom_col()+facet_wrap(~sample)
dev.off()
clonotypes<-as_tibble(ddply(data.trb,"aaSeq",numcolwise(sum)))
clonotypes<-clonotypes[order(-clonotypes$size),c('aaSeq','size')]
topTen<-clonotypes$aaSeq[1:100]
for (j in 1:topTen){
  subset<-data.trb[grepl(topTen[j],data.trb$aaSeq),]
  ddply(subset,"sample",numcolwise(sum))
}
data.trb.topTen<-data.trb[data.trb$aaSeq %in% topTen,]
ggplot(data.trb.topTen,aes(sample,aaSeq))+geom_tile(aes(fill=size),colour = "white")+coord_flip()

