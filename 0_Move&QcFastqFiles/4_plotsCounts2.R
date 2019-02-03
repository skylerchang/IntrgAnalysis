library(tidyverse)
library(here)
library(gsubfn)

targetfolder<-"../../Data/Counts/"
outfolder<-"../../Data/CountPlots/"
outTables<-"../../Data/Tables/"

#parse 'before' names & truncate path
beforeName<-read_table(paste0(targetfolder,"beforeQC-names.txt"),col_names=F)
beforeNameShort<-as.tibble(unlist(strapply(beforeName$X1,"^.*/(.+)$")))
#parse 'before' counts 
beforeCount<-read_table(paste0(targetfolder,"beforeQC-count.txt"),col_names=F)

#parse 'after' names & truncate path
afterName<-read_table(paste0(targetfolder,"afterQC-names.txt"),col_names=F)
afterNameShort<-as.tibble(unlist(strapply(afterName$X1,"^.*/(.+)$")))
#parse 'after' counts 
afterCount<-read_table(paste0(targetfolder,"afterQC-count.txt"),col_names=F)

#add column 'type' with value 'before'; divide count by 4 to account for fastq formoat (4 lines per seq)
type<-rep('before',nrow(beforeNameShort))
before<-as_tibble(cbind(beforeNameShort,beforeCount,type))
colnames(before)<-c('file','count','type')
before$count<-before$count/4

#alternative to previous paragraph
before<-as_tibble(cbind(beforeNameShort,beforeCount))
colnames(before)<-c('file','count')

#add column 'type' with value 'after'; divide count by 4 to account for fastq formoat (4 lines per seq)
type<-rep('after',length(afterNameShort))
after<-as_tibble(cbind(afterNameShort,afterCount,type))
colnames(after)<-c('file','count','type')
after$count<-after$count/4

#alternative to previous paragraph
after<-as_tibble(cbind(afterNameShort,afterCount))
colnames(after)<-c('file','count')

merge(after,before,by='file')

n<-tibble(id=after$file,before=before$count,after=after$count)
n<-n[grep("_R1_",n$id),]
n$id<-sub("_L001_R1_001.fastq.gz","",n$id)
write_excel_csv(n,paste0(outTables,"beforeAfterQC_counts.csv"))

#combine tibbles 'before' & 'after'
all<-as_tibble(rbind(after,before))
all$type<-factor(all$type,levels=c('before','after'))

#generate difference between 'before' & 'after' to use for stacked bargraph
diff<-as_tibble(cbind(before$file,(before$count-after$count),as.character(before$type)))
colnames(diff)<-c('file','count','type')
diff$count<-as.numeric(diff$count)
diff$type<-as.factor(diff$type)

#combine tibbles 'before' & 'after'
all<-as_tibble(rbind(after,diff))
all$type<-factor(all$type,levels=c('before','after'))

dir.create(here(outfolder))

#succeeded by following paragraph
#pdf(here(outfolder,"pairedReads_beforeAfterQC.pdf"), width = 15, height = 25)
#ggplot(all,aes(file,count,fill=type))+geom_col()+coord_flip()+scale_fill_manual(values=c('darkgrey'#,'royalblue3'))
#dev.off()


#limit data set to R1 reads since R1 & R2 are identical
a<-all[grep("R1",all$file),]
a$file<-sub("_L001_R1_001.fastq.gz","",a$file)

pdf(here(outfolder,"pairedReads_beforeAfterQC.pdf"), width = 15, height = 25)
ggplot(a,aes(file,count,fill=type))+geom_col()+coord_flip()+scale_fill_manual(values=c('darkgrey','royalblue3'))
dev.off()



total<-sum(before$count)
