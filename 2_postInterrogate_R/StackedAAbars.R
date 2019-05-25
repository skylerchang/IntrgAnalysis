
library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)
library(plyr)
library(reshape)


######  Read rds file only ######
t<-read_rds('RDS/clntab_vAndJ.rds')
outpath<-'OUT/cdrAaComparison'
loci<-c("TRA","TRB","TRD","TRG","IGH")
n<-100
datalist<-t[[1]]
files_short<-t[[2]]
print(files_short)
######  Choose 2 files from list to be compared ##########
my.sample.1<- readline(prompt = "Enter sample1 name: ")
my.sample.2<- readline(prompt = "Enter sample2 name: ")
sampleID.1<-which(grepl(my.sample.1,t[[2]]))
sampleID.2<-which(grepl(my.sample.2,t[[2]]))
compare_list<-list(datalist[[sampleID.1]],datalist[[sampleID.2]])

plotlist<-list()
for (j in 1:length(loci)){
tablelist<-list()
  for (i in 1:length(compare_list)){
    t<-compare_list[[i]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]
    #subset for locus
    tt<-t[t$vAndJchainSimplified==loci[j],c("aaSeq","size")]
    tt<-tt[!is.na(tt$aaSeq),]
    #collapse lines with identical 'aaSeq' and sum up 'size'
    tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
    tt<-tt[order(-tt$size),]
    #if there are greater than n sequences, collapse sequences of identical length
    if (nrow(tt)>n){
      #convert aaSeqs in string of 'x's
      tt$aaSeq[(n+1):nrow(tt)]<-gsub("[A-Z\\*#]","x",tt$aaSeq[(n+1):nrow(tt)])
      #collapse aaSeq
      tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
      tt$percent<-tt$size/sum(tt$size)
      tt<-tt[order(-tt$size),]
    } else if (nrow(tt)>0 && nrow(tt)<n) {
      tt$percent<-tt$size/sum(tt$size)
      tt<-tt[order(-tt$size),]
    }
 tablelist[[i]]<-tt
  }
####### if there is no empty table list, melt and organize tables into one #######
if ((length(tablelist[[1]])+length(tablelist[[2]])) > 5) {
df1<-as.data.frame(tablelist[[1]])
df1$sample<-c("S1")
df1$size<-NULL
df2<-as.data.frame(tablelist[[2]])
df2$sample<-c("S2")
df2$size<-NULL
df3<-rbind(df1,df2)
df4<-data.frame(df3$aaSeq,df3$sample)
df5<-data.frame(df3$aaSeq,df3$sample,df3$percent)
a<-table(df4)
a1<-as.data.frame.matrix(a)
a2<-cast(df5, df3.aaSeq ~ df3.sample)
dfnew<-data.frame(a2$df3.aaSeq,a2$S1,a2$S2)
dfnew[is.na(dfnew)]<-0
names(dfnew)<-c("aaSeq","S1","S2")
dfnew$length<-nchar(as.character(dfnew$aaSeq))
colourCount <-nrow(dfnew)
# assign color palette to stacked bars
getPalette <- colorRampPalette(brewer.pal(8, "Accent"))
dfnew$color<-getPalette(colourCount)
dfnew$color[grepl("^x.*x$", dfnew$aaSeq)] <-"grey"
dfnew2<-melt(dfnew,id=c("aaSeq","length","color"))
dfnew3<-dfnew2 %>% filter(length < 31)
dfnew3$value<-dfnew3$value*100
# plot
plotname=paste("Aa's frequency in", paste0(loci[j]),"between S1:",paste0(my.sample.1),"and S2:",paste0(my.sample.2))
p<-ggplot(data=dfnew3, aes(y = value, x = variable,fill = aaSeq),borders="white") +  geom_col(position = position_stack(reverse = TRUE)) + 
  scale_fill_manual(values=dfnew3$color) + guides(fill=F) +
  facet_grid( ~ length)+ ggtitle(plotname) + theme_grey(base_size=26)+
  xlab("samples in different Aa lengths") + ylab("percentage %")
plotlist[[j]]<-ggplotGrob(p)
# empty table list will be skipped
}  else  {
   next
}
pdf(paste0(outpath,"/Comparison_",my.sample.1,"VS",my.sample.2,".pdf"),width=26,height=38)
plotlist_new<-compact(plotlist) 
do.call("grid.arrange",c(plotlist_new,ncol=1))
dev.off()
}





