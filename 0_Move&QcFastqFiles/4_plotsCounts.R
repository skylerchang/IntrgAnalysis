library(tidyverse)
library(here)
library(gsubfn)

targetfolder<-"../../Data/Counts/"

beforeName<-read_table(paste0(targetfolder,"before-names.txt"),col_names=F)
beforeNameShort<-as.tibble(unlist(strapply(beforeName$X1,"^.*/(.+)$")))
beforeCount<-read_table(paste0(targetfolder,"before-count.txt"),col_names=F)
afterName<-read_table(paste0(targetfolder,"after-names.txt"),col_names=F)
afterNameShort<-as.tibble(unlist(strapply(afterName$X1,"^.*/(.+)$")))
afterCount<-read_table(paste0(targetfolder,"after-count.txt"),col_names=F)


type<-rep('before',nrow(beforeNameShort))
before<-as_tibble(cbind(beforeNameShort,beforeCount,type))
colnames(before)<-c('file','count','type')
before$count<-before$count/4

type<-rep('after',length(afterNameShort))
after<-as_tibble(cbind(afterNameShort,afterCount,type))
colnames(after)<-c('file','count','type')
after$count<-after$count/4


diff<-as_tibble(cbind(before$file,(before$count-after$count),as.character(before$type)))
colnames(diff)<-c('file','count','type')
diff$count<-as.numeric(diff$count)
diff$type<-as.factor(diff$type)

all<-as_tibble(rbind(after,diff))

ggplot(all,aes(file,count,fill=type))+geom_col()+coord_flip()+scale_fill_manual(values=c('darkgrey','royalblue3'))

total<-sum(before$count)
