library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)
library(plyr)
library(reshape)
library(tidyr)


col1<-brewer.pal(n = 8, name = 'Dark2')
col2<-brewer.pal(n = 12, name = 'Set3')
col3<-brewer.pal(n = 12, name = 'Paired')
col4<-brewer.pal(n = 9, name = 'Pastel1')
col5<-brewer.pal(n = 8, name='Pastel2')
col6<-brewer.pal(n = 8, name = 'Set2')
col7<-brewer.pal(n = 8, name = 'Accent') 
col8<-brewer.pal(n = 9, name = 'Set1')
col9<-brewer.pal(n = 10, name = 'PRGn')
col10<-brewer.pal(n = 10, name = 'RdBu')
col11<-brewer.pal(n = 10, name = 'BrBG')
col<-as.vector(rbind(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11))
col<-unique(col)
col[[78]] <- "#8E0152"
col[[90]] <- "#C51B7D"
col[[81]] <- "#DE77AE"
col[[89]] <- "#F1B6DA"
col[[87]] <- "#3288BD"
col[[82]] <- "#66C2A5"
col[[100]] <-"#D53E4F"
col[[101]] <-"grey40"
col<-unique(col)


d<-read_rds('RDS/clntab_vAndJ.rds')
outpath<-'OUT/CdrAaComparison/'
loci<-c("TRB","IGH")
#Top n clones that are being displayed in separate colors (all other clones are grey)
n<-100
#==========================================================
#datalist contains one tibble per sample
datalist<-d[[1]]
#sample names are stored separately in a vector
files_short<-d[[2]]

for (m in seq(from=1,to=length(datalist),by=2)){
  my.sample.1<- files_short[m]
  my.sample.2<- files_short[m+1]
  name_list<-list(my.sample.1,my.sample.2)
  files_shortA<- gsub("-[AB].*","",my.sample.1)
  full_names<- gsub("_S.","",my.sample.1)
  compare_list<-list(datalist[[m]],datalist[[(m+1)]])
  
  plotlist<-list()
  for (j in 1:length(loci)) {
  tablelist<-list()
  for (i in 1:length(compare_list)){
    t<-compare_list[[i]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]
    #subset for locus
    tt<-t[t$vAndJchainSimplified==loci[j],c("aaSeq","size")]
    #tt<-t[t$vAndJchainSimplified=="TRB",c("aaSeq","size")]
    tt<-tt[!is.na(tt$aaSeq),]
    #collapse lines with identical 'aaSeq' and sum up 'size'
    tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
    tt<-tt[order(-tt$size),]
    tt[tt$aaSeq=="noju",c("aaSeq","size")]
    tt$aaSeq[tt$aaSeq=="noju"]<-"xxxx"
    tt$size[tt$aaSeq=="xxxx"]<-0
    #if there are greater than n sequences, collapse sequences of identical length
    if (nrow(tt)>n){
      #convert aaSeqs in string of 'x's
      tt$aaSeq[(n+1):nrow(tt)]<-gsub("[A-Z\\*#]","x",tt$aaSeq[(n+1):nrow(tt)])
      #collapse aaSeq
      tt$percent<-tt$size/sum(tt$size)
      tt<-tt[order(-tt$size),]
    } else if (nrow(tt)>0 && nrow(tt)<n) {
      tt$percent<-tt$size/sum(tt$size)
      tt<-tt[order(-tt$size),]
    }
    tablelist[[i]]<-tt
  }
  # if there is no empty table list, melt and organize tables into one
  if ((length(tablelist[[1]])+length(tablelist[[2]])) > 5) {
    df1<-as.data.frame(tablelist[[1]])
    if ( nrow(df1) > 100 ) { 
      df1$color<-"grey"
      df1$color[1:100]<-col
    } else {
      df1$color<-"grey"
      df1$color[1:nrow(df1)]<-col[1:nrow(df1)]
    }
    df1$sample<-c("S1")
    df1$size<-NULL
    df1<-as_tibble(ddply(df1,.(aaSeq,color,sample),numcolwise(sum)))
    df1<-df1[order(-df1$percent),]
    
    df2<-as.data.frame(tablelist[[2]])
    if ( nrow(df2) > 100 ) { 
      df2$color<-"grey"
      df2$color[1:100]<-col
    } else {
      df2$color<-"grey"
      df2$color[1:nrow(df2)]<-col[1:nrow(df2)]
    }
    df2$sample<-c("S2") 
    df2$size<-NULL
    df2<-as_tibble(ddply(df2,.(aaSeq,color,sample),numcolwise(sum)))
    df2<-df2[order(-df2$percent),]
    
    df3<-rbind(df1,df2)
    df5<-data.frame(df3$aaSeq,df3$sample,df3$percent)
    df6<-data.frame(df3$aaSeq,df3$sample,df3$color)
    #a2<-cast(df5, df3.aaSeq ~ df3.sample)
    #a3<-cast(df6, df3.aaSeq ~ df3.sample)
    a2<-reshape(df5,direction="wide",idvar="df3.aaSeq",timevar="df3.sample")
    a3<-reshape(df6,direction="wide",idvar="df3.aaSeq",timevar="df3.sample")

    a2<-as.data.frame(a2)
    colnames(a2)<-c("aaSeq","S1","S2")
    a2[is.na(a2)]<-0
    a2$length<-nchar(as.character(a2$aaSeq))
    dfnew2<-melt(a2,id=c("aaSeq","length"))
    
    a3<-as.data.frame(a3)
    colnames(a3)<-c("aaSeq","S1","S2")
    a3[is.na(a3)]<-"grey"
    a3$length<-nchar(as.character(a3$aaSeq))
    dfnew2B<-melt(a3,id=c("aaSeq","length"))
    
    dfnew2<-cbind(dfnew2,subset(dfnew2B,select=value))
    names(dfnew2)<-c("aaSeq","length","variable","value","color")
    dfnew3<-dfnew2 %>% filter(length < 31)
    dfnew3$value<-dfnew3$value*100
    dfnew3$color<-as.character(dfnew3$color)
    
    dfnew4<-subset(dfnew3, value!=0)  
    dfnew4<-dfnew4[order(-dfnew4$value),]
    dfnew4.colored<-dfnew4[!grepl("x",dfnew4$aaSeq),]
    dfnew4.grey<-dfnew4[grepl("x",dfnew4$aaSeq),]
    dfnew4.final<-rbind(dfnew4.colored,dfnew4.grey)
    rownames(dfnew4.final)<-NULL
    dfnew4.final$aaSeq <- reorder(dfnew4.final$aaSeq, dfnew4.final$value)
    dfnew4.final$aaSeq <- factor(dfnew4.final$aaSeq, levels=rev(unique(dfnew4.final$aaSeq)))
    color<-dfnew4.final$color
    names(color) <- dfnew4.final$aaSeq
    # plot
    plotname=paste("Aa's frequency in", paste0(loci[j]),"between S1:",paste0(my.sample.1),"and S2:",paste0(my.sample.2))
    p<-ggplot(data=dfnew4.final, aes(y = value, x = variable, fill = aaSeq),borders="white") +  geom_col(position = position_stack(reverse = FALSE)) + 
      scale_fill_manual(values=color) + guides(fill=F)  + facet_grid( ~ length, switch = "both") + ggtitle(plotname) + theme_grey(base_size=26) + theme(panel.spacing = unit(0, "lines")) + 
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
  if (m > length(datalist)) {
    break
  }   
}











