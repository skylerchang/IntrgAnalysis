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

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ_filtered.rds'
resultsPathDiversity<-'../Results/Diversity/'

#============ functions =======
calculateDiversity<-function(x,y){
  #aggregate all lines with identical aaSeq
  x<-as_tibble(ddply(x[,c("aaSeq","readCount")],"aaSeq",numcolwise(sum)))
  #create tibble with diversity indexes
  tibble(
    filename=y,
    shannon=diversity(x$readCount,index = "shannon",base = exp(1)),
    effectiveSpecies=exp(shannon),
    simpson=diversity(x$readCount,index = "simpson"),
    readCount=sum(x$readCount),
    clonotypeCount=nrow(x)
  )
}

splitFilename<-function(x){
  #transform old control names to required format
  x$filename<-x$filename %>%
    sub("k9-pc-st_D16-07","D16-07-st",.) %>% 
    sub("k9-pc_D16-07","D16-07",.) %>% 
    sub("k9-nt-st_H2O","ntc",.) %>% 
    sub("k9-nt_H2O","ntc-st",.)
  splitFilename<-strsplit(x$filename,"_")
  x$id<-laply(splitFilename, '[[', 1)
  x$ownerPatient<-laply(splitFilename, '[[', 2)
  x$idNumber<-laply(splitFilename, '[[', 3)
  x$pcr<-sub("C[0-9]+P[0-9]+C[0-9]+$","",x$id)
  x$sample<-sub("-D[0-9]+P[0-9]+$","",x$pcr)
  x$submission<-sub("-[0-9a-zA-Z]+$","",x$sample)
  x$replicate<-ifelse(grepl("D[0-9]+P[0-9]*[13579]$",x$pcr),"rep1","rep2")
  if(sampleNoCodesForFraction==T){
    x$fraction<-ifelse(grepl("-1$",sample),"fraction1","fraction2")
  }else{
    x$fraction<-rep("fraction1",length(sample))
  }
  x
}
#===============================

dir.create(resultsPathDiversity,recursive=T)

t<-read_rds(rdsPath)

#datalist contains one tibble per sample
datalist<-t[[1]]
length(datalist)
#sample names are stored separately in a vector
files_short<-t[[2]]
length(files_short)

#initialize tibbles
diversitySummary<-tibble()

for (i in 1:length(files_short)){  
  print(files_short[i])
  a<-datalist[[i]]
  if(is.null(a)){next}
  if(nrow(a)==0){next}
  a$filename<-files_short[i]
  
  #split by locus
  igh<-a[a$locus=='IGH',]
  trb<-a[a$locus=='TRB',]
  
  #calculate diversity indexes and store in 'diversity'
  temp<-calculateDiversity(igh,files_short[i])
  temp$locus<-'IGH'
  diversitySummary<-bind_rows(diversitySummary,temp)
  
  temp<-calculateDiversity(trb,files_short[i])
  temp$locus<-'TRB'
  diversitySummary<-bind_rows(diversitySummary,temp)
}

#===============================================================
#======================= diversity summary =====================
#===============================================================
d<-diversitySummary

#====== split by sample - loci & replicates in rows =============
d.split1<-splitFilename(d)
d.split2<-split(d.split1,d.split1$sample)
#str(d.split2)

#print to xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",d.split1)
for (i in 1:length(d.split2)){
  addWorksheet(wb,d.split2[[i]]$sample[1])
  writeData(wb,d.split2[[i]]$sample[1],d.split2[i],colNames = TRUE)
}
saveWorkbook(wb,paste0(resultsPathDiversity,"diversitySummary.xlsx"),overwrite = T)

#==================== plot - by ownerPatient ===================
#transform d to long format
d.long<-gather(d,index,value,-filename,-locus,-readCount,-clonotypeCount)
d.long<-splitFilename(d.long)

d.long$ownerPatient<-as.factor(d.long$ownerPatient)
colorCount<-nlevels(d.long$ownerPatient)
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))

#jitter all indexes
pdf(paste0(resultsPathDiversity,"diversity_jitter_allIndexes.pdf"))
ggplot(d.long,aes(locus,value))+
  geom_jitter(aes(shape=locus,color=ownerPatient),size=3)+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none")
dev.off()

#by ownerPatient - facet index
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetLocusIndexes.pdf"))
ggplot(d.long,aes(ownerPatient,value))+
  geom_point(aes(shape=locus,color=ownerPatient),size=3)+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#by ownerPatient - facet index - IGH
d.long.igh<-d.long[d.long$locus=='IGH',]
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetIndex_igh.pdf"))
ggplot(d.long.igh,aes(ownerPatient,value))+
  geom_point(aes(shape=locus,color=ownerPatient),size=3)+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1))
dev.off()
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetIndex_igh_box.pdf"))
ggplot(d.long.igh,aes(ownerPatient,value,fill=ownerPatient))+
  geom_boxplot()+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#by ownerPatient - facet index - TRB
d.long.trb<-d.long[d.long$locus=='TRB',]
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetIndex_trb.pdf"))
ggplot(d.long.trb,aes(ownerPatient,value))+
  geom_point(aes(shape=locus,color=ownerPatient),size=3)+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1))
dev.off()
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetIndex_trb_box.pdf"))
ggplot(d.long.trb,aes(ownerPatient,value,fill=ownerPatient))+
  geom_boxplot()+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#==================== plot - by readCount ===================
#point - by readCount - all indexes 
pdf(paste0(resultsPathDiversity,"readCountVsDiversity_facetLocusIndex.pdf"))
ggplot(d.long,aes(readCount,value))+
  geom_point(mapping=aes(shape=locus,color=ownerPatient),size=3)+
  facet_grid(index~locus,scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none")
dev.off()

#========== transpose locus to column - replicates by line =============
d2<-d %>%
  gather(key, value, -filename, -locus) %>%
  unite(indexByLocus, key, locus) %>%
  spread(indexByLocus,value)

d2<-splitFilename(d2)

#save as rds file
saveRDS(d2,paste0(resultsPathDiversity,"diversitySummary.rds"))

#================= merge replicates ================
d2<-splitFilename(d)
d2<-subset(d2,select=c(submission,sample,locus,replicate,shannon,effectiveSpecies,simpson))
#d2<-subset(d2,select=-c(filename,idNumber,id,pcr))
colnames(d2)

d2<-d2 %>% gather(index,value,shannon:simpson) %>%
  unite(index2,index,replicate) %>%
  spread(index2,value)
subset(d2,select=c(fraction,replicate))

d2$shannon.mean<-(d2$shannon_rep1+d2$shannon_rep2)/2
d2$simpson.mean<-(d2$simpson_rep1+d2$simpson_rep2)/2
d2$effectiveSpecies.mean<-(d2$effectiveSpecies_rep1+d2$effectiveSpecies_rep2)/2

d2<-subset(d2,select=c(submission,sample,locus,shannon.mean,simpson.mean,effectiveSpecies.mean))

d2<-d2 %>% gather(index,value,-submission,-sample,-locus) %>%
  unite(index3,index,locus) %>%
  spread(index3,value)

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,'summary')
writeData(wb,'summary',d2)
saveWorkbook(wb,paste0(resultsPathDiversity,'diversity_withMeans.xlsx'),overwrite = T)

#save as rds
write_rds(d2,paste0(resultsPathDiversity,'diversity_withMeans.rds'))
