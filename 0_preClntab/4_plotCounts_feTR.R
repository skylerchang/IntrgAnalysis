library(tidyverse)
library(here)
library(gsubfn)
library(reshape2)

targetfolder<-"../../Data/Counts/"
outfolder<-"../../Data/CountPlots/"
outTables<-"../../Data/Tables/"

dir.create(outfolder)

#parse 'before' names & truncate path
name<-read_table(paste0(targetfolder,"beforeQC-names.txt"),col_names=F)
#parse 'before' counts 
beforeCount<-read_table(paste0(targetfolder,"beforeQC-count.txt"),col_names=F)
#parse 'after' counts 
afterCount<-read_table(paste0(targetfolder,"afterQC-count.txt"),col_names=F)

#bind columns, removed duplicate lines (R1/R2), truncate file name
d<-bind_cols(name,beforeCount,afterCount)
colnames(d)<-c('file','before','after')
d<-d[grep("_R1_",d$file),]
d$file<-sub("_L001_R1_001.fastq","",d$file)
d$file<-basename(d$file)


#melt
dd<-melt(d,id='file')
dd$tissue<-NA
dd$tissue[grepl("L_",dd$file)]<-'Lymph node'
dd$tissue[grepl("S_",dd$file)]<-'Spleen'
dd$tissue[grepl("T_",dd$file)]<-'Thymus'
colnames(dd)<-c('file','type','count','tissue')

dd<-as_tibble(dd)

#bargraph before/after - by file
pdf(paste0(outfolder,"preClntab_beforeAfterQc_bar_byFile.pdf"))
ggplot(dd,aes(file,count,fill=type))+geom_bar(position="dodge", stat="identity")+coord_flip()+scale_fill_brewer(palette="Accent")
dev.off()

#boxplot before/after - by tissue
pdf(paste0(outfolder,"preClntab_beforeAfterQC_box_bySample.pdf"))
ggplot(dd,aes(tissue,count,fill=type))+geom_boxplot()+coord_flip()+geom_point(pch = 21,position = position_jitterdodge(jitter.width = 0.3))+scale_fill_brewer(palette='Dark2')
dev.off()




