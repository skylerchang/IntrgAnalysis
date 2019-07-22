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

source("2_postInterrogate_R/functions.R")

#=========== adjust =============
sampleNoCodesForFraction<-F        #T or F
run<-28

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

#==================== plot - by ownerPatient ===================
#transform d to long format
d.long<-gather(d,index,value,-filename,-locus,-readCount,-clonotypeCount)
d.long<-splitFilename(d.long,sampleNoCodesForFraction,run)

d.long$ownerPatient<-as.factor(d.long$ownerPatient)
colorCount<-nlevels(d.long$ownerPatient)
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))

#jitter all indexes
pdf(paste0(resultsPathDiversity,"diversity_jitter_allIndexes.pdf"))
ggplot(d.long,aes(locus,value))+
  geom_jitter(aes(shape=fraction,color=ownerPatient),size=3)+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="none")
dev.off()

d.long$fraction
#by ownerPatient - facet index
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetLocusIndexes.pdf"))
ggplot(d.long,aes(ownerPatient,value))+
  geom_point(aes(shape=fraction,color=ownerPatient),size=3)+
  facet_grid(index~locus, scales="free")+
  scale_color_manual(values=getPalette(colorCount))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1))
dev.off()

#by ownerPatient - facet index - IGH
d.long.igh<-d.long[d.long$locus=='IGH',]
pdf(paste0(resultsPathDiversity,"diversity_byOwnerPatient_facetIndex_igh.pdf"))
ggplot(d.long.igh,aes(ownerPatient,value))+
  geom_point(aes(shape=replicate,color=ownerPatient),size=3)+
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

d2<-splitFilename(d2,sampleNoCodesForFraction,run)

#save as xlsx file
wb<-createWorkbook()
addWorksheet(wb,"all")
writeData(wb,"all",d2)
saveWorkbook(wb,paste0(resultsPathDiversity,"diversity_locusAsColumn.xlsx"),overwrite = T)

#save as rds file
saveRDS(d2,paste0(resultsPathDiversity,"diversity_locusAsColumn.rds"))

#============= transpose reps to column - locus by line -> merge replicates ================
d2<-splitFilename(d,sampleNoCodesForFraction,run)
d2<-subset(d2,select=c(submission,locus,replicate,fraction,readCount,clonotypeCount,shannon,effectiveSpecies,simpson))
colnames(d2)

#omit controls for now (quick & dirty fix) -> go back and give unique submission name
d2<-d2[!grepl("k9",d2$submission),]

d2<-d2 %>% gather(index,value,readCount:simpson) %>%
  unite(index2,index,replicate) %>%
  spread(index2,value)

d2$readCount.mean<-(d2$readCount_rep1+d2$readCount_rep2)/2
d2$clonotypeCount.mean<-(d2$clonotypeCount_rep1+d2$clonotypeCount_rep2)/2
d2$shannon.mean<-(d2$shannon_rep1+d2$shannon_rep2)/2
d2$simpson.mean<-(d2$simpson_rep1+d2$simpson_rep2)/2
d2$effectiveSpecies.mean<-(d2$effectiveSpecies_rep1+d2$effectiveSpecies_rep2)/2

d2<-subset(d2,select=-c(readCount_rep1,readCount_rep2,clonotypeCount_rep1,clonotypeCount_rep2,shannon_rep1,shannon_rep2,simpson_rep1,simpson_rep2,effectiveSpecies_rep1,effectiveSpecies_rep2))

d2<-d2 %>% gather(key,value,-submission,-locus,-fraction) %>%
  unite(index3,key,locus) %>%
  spread(index3,value)

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,'all')
writeData(wb,'all',d2)
saveWorkbook(wb,paste0(resultsPathDiversity,'diversity_withMeans.xlsx'),overwrite = T)

#save as rds
write_rds(d2,paste0(resultsPathDiversity,'diversity_withMeans.rds'))
