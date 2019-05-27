### ComputeCanada Version needs to load modules####
module load r/3.5.0
module load gcc
module load openmpi

### Install R packages before running script ######
#install.packages("tidyverse")
#install.packages("here")
#install.packages("gsubfn")
#install.packages("reshape2")
#install.packages("proto")




library(tidyverse)
library(here)
library(gsubfn)
library(reshape2)
#library(proto)


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
d$file<-sub("_L001_R1_001.fastq.gz","",d$file)
d$file<-basename(d$file)


#melt
dd<-melt(d,id='file')
dd$rep<-ifelse(grepl("-[0-9]+[A-Z]*D[0-9]+P[0-9]*[13579]C",dd$file),'rep1','rep2')

#dd$sample<-sub('[0-9]D[0-9].*P[0-9].*C[0-9].*P[0-9].*C[0-9]_', '', dd$file)
dd$sample<-sub("D[0-9]+P[0-9]+C[0-9]+P[0-9]+C[0-9]+$","",dd$file)
dd$fraction<-ifelse(grepl("1$",dd$sample),'pellet','supernatant')
#dd$submission<-sub('-[0-9].*-[A-Z].*[a-z]', '', dd$sample)
dd$submission<-sub("-[0-9]+[A-Z]*$","",dd$sample)
colnames(dd)<-c('file','type','count','rep','sample','fraction','submission')
dd<-as_tibble(dd)

#bargraph before/after - by file
pdf(paste0(outfolder,"preClntab_beforeAfterQc_bar_byFile.pdf"))
ggplot(dd,aes(file,count,fill=type))+geom_bar(position="dodge", stat="identity")+coord_flip()
dev.off()

#boxplot before/after - by sample
pdf(paste0(outfolder,"preClntab_beforeAfterQC_box_bySample.pdf"))
ggplot(dd,aes(sample,count,fill=type))+geom_boxplot()+coord_flip()
dev.off()

#boxplot before/after - by sample - log x axis
pdf(paste0(outfolder,"preClntab_beforeAfterQC_box_bySample_log.pdf"))
ggplot(dd,aes(sample,log(count),fill=type))+geom_boxplot()+coord_flip()
dev.off()

pdf(paste0(outfolder,"preClntab_countByFraction.pdf"))
ggplot(dd,aes(submission,count,fill=fraction))+geom_boxplot()+coord_flip()
dev.off()

pdf(paste0(outfolder,"preClntab_countByFraction_log.pdf"))
ggplot(dd,aes(submission,log(count),fill=fraction))+geom_boxplot()+coord_flip()
dev.off()


