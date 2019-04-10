#this script is specific to Run02; save under new name and adapt to other runs)


library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)

#========== adjust the following variables ================
t<-read_rds('../../Data/Clntab_RDS/clntab_vAndJ.rds')
outpath<-'../../Results/Plots/cdrAaLengthVsAaSeq'
loci<-c("TRA","TRB","TRD","TRG","IGH")
#Top n clones that are being displayed in separate colors (all other clones are grey)
n<-80

#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]


for (i in 1:length(datalist)){
  t<-datalist[[i]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]

  nrow(t)
  #================ aaLength vs. aaSeq =================================
  #collapse reads with identical aaSeq
  library(plyr)
  
  plotlist<-list()
  for (j in 1:length(loci)){
    #subset for locus
    tt<-t[t$vAndJchainSimplified==loci[j],c("aaSeq","size")]
    tt<-tt[!is.na(tt$aaSeq),]
    #collapse lines with identical 'aaSeq' and sum up 'size'
    tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
    #sort in descending order of 'size'
    tt<-tt[order(-tt$size),]
    #assign color - create color vector
    #b<-rep(brewer.pal(8,"Accent"),ceiling(nrow(tt)/8))
    #tt$color<-b[1:nrow(tt)]
    
    #if there are greater than n sequences, collapse sequences of identical length
    if (nrow(tt)>n){
      #convert aaSeqs in string of 'x's
      tt$aaSeq[(n+1):nrow(tt)]<-gsub("[A-Z\\*#]","x",tt$aaSeq[(n+1):nrow(tt)])
      #collapse aaSeq
      tt<-as_tibble(ddply(tt,"aaSeq",numcolwise(sum)))
      #assign color
      b<-rep(brewer.pal(8,"Accent"),ceiling(nrow(tt)/8))
      tt$color<-b[1:nrow(tt)]
      tt$color[(n+1):nrow(tt)]<-"grey"
    } else if (nrow(tt)>0 && nrow(tt)<n) {
            b<-rep(brewer.pal(8,"Accent"),ceiling(nrow(tt)/8))
            tt$color<-b[1:nrow(tt)]
    } else if (nrow(tt)==0) {
            next
    }
    
    #create length column
    tt$length<-nchar(tt$aaSeq)
    #plot
    plotname=paste(paste0(loci[j]),"Aa's length vs. size")
    p<-ggplot(tt,aes(length,size,fill=aaSeq))+geom_col(position = position_stack(reverse = TRUE))+guides(fill=F)+scale_fill_manual(values=tt$color)+xlim(1,35)+ggtitle(plotname)
    plotlist[[j]]<-ggplotGrob(p)
  }
  pdf(paste0(outpath,"/cdrAaSeqByAaLength_",files_short[i],"_size25_unicolor.pdf"))
  plotlist_new<-compact(plotlist) 
  do.call("grid.arrange",c(plotlist_new,ncol=2,top=files_short[i]))
  #marrangeGrob(plotlist,nrow=2,ncol=2)
  dev.off()
}

