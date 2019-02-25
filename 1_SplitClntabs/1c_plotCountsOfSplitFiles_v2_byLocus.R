library(tidyverse)
library(RColorBrewer)
library(here)


dir<-"../../Results/Clntab/"

#****** 'uniread' file has leading zeros -> adjust script or remove manually (import fails)****

for (type in c('reads','unireads')){
  t<-read_delim(paste0(dir,'clntab-',type,'.txt'), delim = " ", col_names = F)
  colnames(t)<-c('count','filename')
  
  #generate id from file name
  t$id<-basename(t$filename)
  t$id<-sub("_L001_R1_001.*", "",t$id)
  
  t$organ<-rep(NA,nrow(t))
  t$organ[grepl('L',t$filename)]<-'LymphNode'
  t$organ[grepl('T',t$filename)]<-'Thymus'
  t$organ[grepl('S',t$filename)]<-'Spleen'
  t$organ<-as.factor(t$organ)
  
  pdf(paste0('Plots/',type,'.pdf'))
  ggplot(t,aes(id,count,fill=organ))+geom_col()+coord_flip()+scale_fill_brewer(palette = "Dark2")
  dev.off()
}

