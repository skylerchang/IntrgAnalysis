#-one sample per file (i.e. duplicates are in separate files); 
#-multiple loci per file
#-no template filtering



library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)

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
t<-read_rds('../Data/Clntab_RDS/clntab_vAndJ.rds')
outpath<-'../Results/CdrAaLengthVsAaSeq/'
loci<-c("TRA","TRB","TRD","TRG","IGH")
#Top n clones that are being displayed in separate colors (all other clones are grey)
n<-100

#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

dir.create(outpath<-'../Results/CdrAaLengthVsAaSeq/',recursive=T)

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
        #if there are greater than n sequences, collapse sequences of identical length
        if (nrow(tt)>n){
            #convert aaSeqs in string of 'x's
            tt$aaSeq[(n+1):nrow(tt)]<-gsub("[A-Z\\*#]","x",tt$aaSeq[(n+1):nrow(tt)])
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
        p<-ggplot(tt,aes(length,size,fill=aaSeq))+geom_col(position = position_stack(reverse = FALSE))+guides(fill=F)+scale_fill_manual(values=tt$color)+xlim(3,30)+ggtitle(plotname)
        plotlist[[j]]<-ggplotGrob(p)
    }
    pdf(paste0(outpath,"/cdrAaSeqByAaLength_",files_short[i],"_size25_unicolor.pdf"))
    plotlist_new<-compact(plotlist) 
    do.call("grid.arrange",c(plotlist_new,ncol=2,top=files_short[i]))
    dev.off()
}
