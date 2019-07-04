library(tidyverse)
library(plyr)
library(here)

setwd(here())
getwd()

#========== adjust variables ================
infile<-'../Data/clntab_RDS/clntab_vAndJ.rds'
outpath<-'../Results/RepertoireOverlap/'

#==========================================================
t<-read_rds(infile)
dir.create(outpath,recursive=T)

#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

pcrId<-sub("_.*","",files_short)
sampleId<-sub("-D[0-9]+P[0-9]+$","",pcrId)

data<-tibble()

for (i in 1:length(datalist)){
  t<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]
  t$pcrId<-pcrId[i]
  t$sampleId<-sampleId[i]
  data<-bind_rows(data,t)
}

# filter by "^C.{3,30}[FW]$" & locus
d<-data[grepl("^C.{3,30}[FW]$",data$aaSeq),]
d<-d[!d$vAndJchainSimplified=='TRA',]
#create col replicate
d$replicate<-ifelse(grepl("[13579]$",d$pcrId),"rep1","rep2")

#split by locus
igh<-d[d$vAndJchainSimplified=='IGH',]
trb<-d[d$vAndJchainSimplified=='TRB',]

#print size distribution for igh & trb individually
ggplot(d,aes(aaLength,size))+geom_col()+facet_wrap(~vAndJchainSimplified)

#============ aalength by sample - igh =================
pdf(paste0(outpath,'facet-igh.pdf'))
ggplot(igh,aes(aaLength,size))+geom_col()+facet_wrap(~sample)
dev.off()

#============ aaLength by sample - trb =================
pdf(paste0(outpath,'facet-trb.pdf'))
ggplot(trb,aes(aaLength,size))+geom_col()+facet_wrap(~sample)
dev.off()

#++++++++++++++++++ function - determine overlap ++++++++++++
determineOverlap<-function(x){
  samples<-unique(x$sampleId)
  overlap<-tibble()
  
  for (i in 1:length(samples)){
    t<-x[x$sampleId==samples[i],c("aaSeq","size","replicate")]
    t<-as_tibble(ddply(t,c("aaSeq","replicate"),numcolwise(sum)))
    t<-t[t$size>=10,]
    freq<-as_tibble(count(t,vars=c("aaSeq")))
    o<-freq$aaSeq[freq$freq>1]
    temp<-tibble(clone=o,sampleId=samples[i])
    overlap<-bind_rows(overlap,temp)
  }
  
  f<-as_tibble(count(overlap,vars=c("clone")))
  f[f$freq>1,]
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

summary.igh<-determineOverlap(trb)
summary.trb<-determineOverlap(igh)

capture.output(summary.igh,file=paste0(outpath,'summary_igh.txt'))
capture.output(summary.trb,file=paste0(outpath,'summary_trb.txt'))


#====== this is (non-serious) run specific code (histiocytoma), variables are no longer up to date ======================

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

