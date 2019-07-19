#script calculates diversity for split 'Clntab_RDS' files: 'Clntab_RDS_noTemplate' & 'Clntab_RDS_templateOnly'

library(tidyverse)
library(plyr)
library(here)
library(reshape2)
library(openxlsx)
library(RColorBrewer)
library(ggforce)

setwd(here())
getwd()

#=========== adjust =============
loci<-c("IGH","TRB")
run<-28
sampleNoCodesForFraction<-F        #T or F

rdsPath<-'../Data/Clntab_RDS/clntab_vAndJ_filtered.rds'
resultsPathLocus<-'../Results/LocusStats/'

#===============================
dir.create(resultsPathLocus,recursive=T)

t<-read_rds(rdsPath)

#datalist contains one tibble per sample
datalist<-t[[1]]
length(datalist)
#sample names are stored separately in a vector
files_short<-t[[2]]
length(files_short)

#initialize things
locusSummary<-tibble()

for (i in 1:length(files_short)){  
  print(files_short[i])
  a<-datalist[[i]]
  if(is.null(a)){
    print("no data - NULL")
    next
  }
  count.all<-nrow(a)
  if(count.all==0){
    print("no data")
    next
  }
  a$filename<-files_short[i]
  
  #============ create locus summary =======
  igh<-a[a$locus=='IGH',]
  trb<-a[a$locus=='TRB',]
  
  temp<-tibble(
    filename=files_short[i],
    totalReads=sum(a$readCount),
    igh.count=sum(igh$readCount),
    igh.percent=igh.count/totalReads*100,
    trb.count=sum(trb$readCount),
    trb.percent=trb.count/totalReads*100,
    locusRatio=igh.count/trb.count
    )
  locusSummary<-bind_rows(locusSummary,temp)
}

#===============================================================
#========================== locus summary ======================
#===============================================================
l<-locusSummary

#create additional columns
l<-splitFilename(l,sampleNoCodesForFraction,run)

#=============== save as xlsx & rds =================
#save as rds
saveRDS(l,paste0(resultsPathLocus,"locusSummary.rds"))

#save as xlsx
wb<-createWorkbook()
addWorksheet(wb,"summary")
writeData(wb,"summary",l)
saveWorkbook(wb,paste0(resultsPathLocus,"locusSummary.xlsx"),overwrite = T)

#=============== plot locus count data =================
#subset count data, reshape to long format, split file name
l.count.long<-subset(l,select=c(filename,igh.count,trb.count)) %>%
  gather(locus,value,-filename) %>%
  splitFilename(.,sampleNoCodesForFraction,run)
#create 4 replicates by incorporating the locus information -> plot loci in different colors
l.count.long$replicate2<-paste(l.count.long$locus,l.count.long$replicate,sep='-')

pdf(paste0(resultsPathLocus,"ighVsTrb_withReplicates_barplot.pdf"))
ggplot(l.count.long,aes(locus,value,fill=replicate2))+geom_col(position=position_dodge())+facet_wrap(~sample)+scale_fill_brewer(palette="Paired")
dev.off()

#================= plot locus ratio data =================
l.splitFilename<-splitFilename(l,sampleNoCodesForFraction,run)

#total reads vs locus ratio
pdf(paste0(resultsPathLocus,"readsVsLocusRatio_all.pdf"))
ggplot(l.splitFilename,aes(totalReads,locusRatio,fill=sample))+geom_point(pch=21)+theme(legend.position="top")
dev.off()

#total reads vs locus ratio - with label
pdf(paste0(resultsPathLocus,"readsVsLocusRatio_all_withLabel.pdf"))
ggplot(l.splitFilename,aes(totalReads,log(locusRatio),fill=sample))+geom_point(pch=21)+geom_label(aes(label=sample))+theme(legend.position="none")
dev.off()

pdf(paste0(resultsPathLocus,"readsVsSample_all_boxplot.pdf"))
ggplot(l.splitFilename,aes(totalReads,log(locusRatio),group=sample))+geom_boxplot(aes(fill=sample))
dev.off()

pdf(paste0(resultsPathLocus,"readsVslocusRatio_all_smoothedCondMeans.pdf"))
ggplot(l.splitFilename,aes(totalReads,log(locusRatio)))+geom_point()+geom_smooth()
dev.off()

#remove outliers 
#exclude<-c("k9-nt-H2O-0","k9-nt-st-H2O-1","09-020952-16")
#l.splitFilename<-l.splitFilename[!l.splitFilename$sample %in% exclude,]

