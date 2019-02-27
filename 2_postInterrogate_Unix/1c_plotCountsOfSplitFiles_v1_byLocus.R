library(tidyverse)
library(RColorBrewer)


setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Normal/Bioinformatics/Analysis_Nikos/Clntab_2015/Results')

#****** 'uniread' file has leading zeros -> adjust script or remove manually (import fails)****

for (type in c('reads','unireads')){
  t<-read_delim(paste0('clntab-',type,'.txt'), delim = " ", col_names = F)
  colnames(t)<-c('count','filename')
  
  #shorten file name
  t$id<-sub(".*Data/", "",t$filename)
  t$id<-sub("/Cat.*", "",t$id)
  
  t$organ<-rep(NA,nrow(t))
  t$organ[grepl('Ln',t$filename)]<-'LymphNode'
  t$organ[grepl('Th',t$filename)]<-'Thymus'
  t$organ[grepl('Sp',t$filename)]<-'Spleen'
  t$organ<-as.factor(t$organ)
  
  pdf(paste0('Plots/',type,'.pdf'))
  ggplot(t,aes(id,count,fill=organ))+geom_col()+coord_flip()+scale_fill_brewer(palette = "Dark2")
  dev.off()
}

