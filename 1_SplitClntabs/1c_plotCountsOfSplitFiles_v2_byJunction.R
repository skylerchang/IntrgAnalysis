#takes clntab files that have been split into junction and noju; no prior separation by locus


library(tidyverse)
library(RColorBrewer)
library(here)



path<-"../../Data/Clntab/"

dirs<-list.files(path)

for (i in 1:length(dirs)){
  f<-list.files(paste0(path,dirs[i]),pattern = "junc.clntab")
  t<-read_delim(paste0(path,dirs[i],"/",f), delim = "\t", col_names = T)
  colnames(t)<-c('five','three','aaSeq','size')
  
  t$locus<-NA
  t$locus[grepl("IGH",t$three)]<-'igh'
  t$locus[grepl("IGK",t$three)]<-'igk'
  t$locus[grepl("TRA",t$three)]<-'tra'
  t$locus[grepl("TRB",t$three)]<-'trb'
  t$locus[grepl("TRD",t$three)]<-'trd'
  t$locus[grepl("TRG",t$three)]<-'trg'
  t$locus<-as.factor(t$locus)
  
  ggplot(t,aes(locus))+geom_bar()
}

