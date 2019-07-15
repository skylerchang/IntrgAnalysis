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
run<-30
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
  count.all<-nrow(a)
  if(nrow(a)==0){next}
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

#function splits file name into components creating a new col for each; starting format:
#<id>_<ownerPatient>_<idNumber>
#strsplit by '_' -> 3 parts: 1) id, 2) ownerPatient, 3) idNumber
#1) id is then broken down further into: pcr, sample, submission, replicate, fraction
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

#create additional columns
l<-splitFilename(l)

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
  splitFilename()
#create 4 replicates by incorporating the locus information -> plot loci in different colors
l.count.long$replicate2<-paste(l.count.long$locus,l.count.long$replicate,sep='-')

pdf(paste0(resultsPathLocus,"ighVsTrb_withReplicates_barplot.pdf"))
ggplot(l.count.long,aes(locus,value,fill=replicate2))+geom_col(position=position_dodge())+facet_wrap(~sample)+scale_fill_brewer(palette="Paired")
dev.off()

#================= plot locus ratio data =================
l.splitFilename<-splitFilename(l)

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

