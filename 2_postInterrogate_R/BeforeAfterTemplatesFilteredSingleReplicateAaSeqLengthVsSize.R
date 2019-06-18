library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)
library(plyr)


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
col<-c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11)
col<-unique(col)


setwd(here())
getwd()

#========== adjust the following variables ================
t<-read_rds('../Clntab_RDS/clntab_vAndJ.rds')
outpath<-'../Results/CdrAaLengthTemplate/'
loci<-c("IGH")
#Top n clones that are being displayed in separate colors (all other clones are grey)
n<-100

#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

dir.create(outpath<-'../Results/CdrAaLengthTemplate/',recursive=T)

for (i in 1:length(datalist)){
  a <- datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","completeNtSeq","vAndJchainSimplified")]
  a <- subset(a,(a$aaSeq!="noju"))
  a.all <- a
  #============ check for templates =======
  a$template<-'no'
  #Akash‘s templates
  a$template[grepl("TGTGCATCACGACACAGTGGTCTGG",a$completeNtSeq)]<-'t1'
  a$template[grepl("TGTGCATCACGACCAGATCCACAGATCCATTGGTTACTGG",a$completeNtSeq)]<-'t2'
  #Tamara‘s templates
  a$template[grepl("TGTTCGCCTTATCGCCTTATGG",a$completeNtSeq)]<-'IGHV1-30_IGHJ4'
  a$template[grepl("TGTCTAGTACGCCTCTCTGCCTCTCTGCTAGTACGTGG",a$completeNtSeq)]<-'IGHV1-30_IGHJ6'
  a$template[grepl("TGTGCTTCTGCCTTTCTGCCTGCTCAGGATTCTGCCTGCTCAGGAGCTCAGGATTCTCTGG",a$completeNtSeq)]<-'IGHV3-1_IGHJ4'
  a$template[grepl("TGTGAGGAGTCCGTAGAGAGAGGAGTCCAGCGTAGCCATGCCTAAGGAGTCCCAGCCTCGGTAGAGAGAGCGCTGG",a$completeNtSeq)]<-'IGHV3-1_IGHJ6'
  a$template[grepl("TGTAAGGAGTAACTGCATAACTGCATACTAAGCCTAAGGAGTATGG",a$completeNtSeq)]<-'IGHV4-1_IGHJ4'
  a$template[grepl("TGTGTGCCTCTTTCCTCTACTAGATCGCCTCTCTATTATCCTCTAGAGTAGAGTAAGGAGTAGATCGCTATCCTCTGTAAGGAGTCCTCTACCTGG",a$completeNtSeq)]<-'IGHV4-1_IGHJ6'
  a.woTemp<-a[a$template=='no',]
  a.tempOnly<-a[a$template!='no',]
  
  #============ summarize templates for sample and store in tibble =======
  temp<-as_tibble(ddply(a.tempOnly[,c("template","size")],"template",numcolwise(sum)))
  if(nrow(temp)==0){next}
  temp$sample<-files_short[i]
  templateSummary<-bind_rows(templateSummary,temp)
 ####### make compared Plots between all reads AAseqLength VS Size and without template reads ######  
  plotlist<-list()
  for (j in 1:length(loci)){
    create.stackedbar<-function(t) {
    #subset for locus
    tt<-t[t$vAndJchainSimplified==loci[j],c("aaSeq","size")]
    tt<-tt[!is.na(tt$aaSeq),]
    #collapse lines with identical 'aaSeq' and sum up 'size'
    tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
    #sort in descending order of 'size'
    tt<-tt[order(-tt$size),]
    #if there are greater than n sequences, collapse sequences of identical length
    if (nrow(tt)>n){
      #convert aaSeqs in string of 'x's
      tt$aaSeq[(n+1):nrow(tt)]<-gsub("[A-Z\\*#]","x",tt$aaSeq[(n+1):nrow(tt)])
      #collapse aaSeq
      #assign color
      tt$color<-"grey"
      tt$color[1:100]<-col
      tt<-as_tibble(ddply(tt,.(aaSeq,color),numcolwise(sum)))
    } else if (nrow(tt)>0 && nrow(tt)<n) {
      tt$color<-"grey"
      tt$color[1:nrow(tt)]<-col[1:nrow(tt)]
    } else if (nrow(tt)==0) {
      next
    }
    #create length column
    tt$length<-nchar(tt$aaSeq)
    tt$aaSeq <- reorder(tt$aaSeq, tt$size)
    tt$aaSeq <- factor(tt$aaSeq, levels=rev(unique(tt$aaSeq)))
    colors<-tt$color
    names(colors) <- tt$aaSeq
    #plot
    plotname=paste(paste0(loci[j]),"Aa's length vs. size")
    p<-ggplot(tt,aes(length,size,fill=aaSeq))+geom_col(position = position_stack(reverse = FALSE))+guides(fill=F)+scale_fill_manual(values=colors)+xlim(1,35)+ggtitle(plotname)
    #plotlist[[j]]<-ggplotGrob(p)
    }
    plot1<-create.stackedbar(a.all)
    plot2<-create.stackedbar(a.woTemp)
    plotvs<-list(plot1,plot2)
    for (m in 1:length(plotvs)) {
      plotlist[[m]]<-ggplotGrob(plotvs[[m]])  
    }
  }
  pdf(paste0(outpath,"/TemplatesVsRemoved_",files_short[i],".pdf"))
  plotlist_new<-compact(plotlist) 
  do.call("grid.arrange",c(plotlist_new,ncol=2,top=paste(paste0(files_short[i],": Templates(w) vs Removed(wo)"))))
  dev.off()
}
