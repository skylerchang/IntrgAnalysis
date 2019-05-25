library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(seqinr)
library(plyr)
library(igraph)
library(stringdist)


### run color100 first ######
col1<-brewer.pal(n = 8, name = 'Dark2')
col2<-brewer.pal(n = 12, name = 'Set3')
col3<-brewer.pal(n = 12, name = 'Paired')
col4<-brewer.pal(n = 9, name = 'Pastel1')
col5<-brewer.pal(n = 8, name='Pastel2')
col6<-brewer.pal(n = 8, name = 'Set2')
col7<-brewer.pal(n = 8, name = 'Accent') 
col8<-brewer.pal(n = 9, name = 'Set1')
col9<-brewer.pal(n = 11, name = 'PRGn')
col10<-brewer.pal(n = 11, name = 'RdBu')
col11<-brewer.pal(n = 11, name = 'BrBG')
col<-as.vector(rbind(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11))
col<-unique(col)
color100<- col[1:100]

####  this version is to read clntab files only #####
targetFolder<-"../../Clntab/"
outFolder<-"../../Result/"
#### e.g. Run25 S1 and S2 ######
filesA<-list.files(targetFolder,pattern ="16-051749",recursive = T)

loci<-c("B","H")
aanum<-100
for (i in 1:length(loci)) {
  pdf(paste0(outFolder,"/cln2igraph_",loci[i],"_loci.pdf"),width=9,height=11)
  par(mfrow=c(4,4),mar=c(1,1,1,1))
  for (j in seq(from=8,to=23,by=1))  {
    for (k in 1:length(filesA)) {
     t<-read_tsv(paste0(targetFolder,filesA[1]))
      files_shortA <- basename(sub("_L001_R1_001.fastq.processed.junctioned.profiled.clntab","",filesA[k]))
      ##### AB only for 16-051749 clntab files #####
      files_shortA <- gsub(".*[AB]_","",files_shortA)
      d<-subset(t, select=c("sequence.JUNCTION.nt seq","sequence.JUNCTION.nt seq.len","sequence.JUNCTION.aa seq","sequence.JUNCTION.aa seq.len","sequence.size","sequence.chain","sequence.functionality"))
      colnames(d)<-c("DNASeq","DNALength","aaSeq","aaLength","size","chain","functions")
      dd<-d[d$chain==loci[i],c("aaSeq","size","aaLength")]
      tt2<-as_tibble(ddply(dd,.(aaSeq,aaLength),numcolwise(sum)))
      tt2<-tt2[order(-tt2$size),]
      tt2<-tt2[complete.cases(tt2),]
      tt2$color<-"grey" 
      ###### assign 100 color ######  
      tt2$color[1:100]<-color100
      test<-tt2[tt2$aaLength==j,c("aaSeq","size","color")]
      test<-as.data.frame(test)
      test<-test[complete.cases(test),]
      test$percent<-test$size/sum(tt2$size)
      if ( nrow(test) > aanum ) {
        test2<-head(test,n=aanum)
        uniqueids2 <- as.character(test2$aaSeq)
        hdids2 <- stringdistmatrix(uniqueids2,uniqueids2,method = "hamming")  
        rownames(hdids2) <- uniqueids2
        colnames(hdids2) <- uniqueids2
        links2 <- hdids2
        net2 <- graph_from_adjacency_matrix(links2<3)
        net2 <- simplify(net2, remove.loops = T)
        net2 <- as.undirected(net2)
        plot(net2,layout=layout.fruchterman.reingold(net2),vertex.label=NA,vertex.color=test2$color,vertex.size=test2$percent*10000)
        title(main=paste0(files_shortA,"len:",j,"Aas"))
      } else if (nrow(test) > 0 && nrow(test) < aanum) { 
        uniqueids <- as.character(test$aaSeq)
        hdids <- stringdistmatrix(uniqueids,uniqueids,method = "hamming")  
        rownames(hdids) <- uniqueids
        colnames(hdids) <- uniqueids
        links <- hdids
        net1 <- graph_from_adjacency_matrix(links<3)
        net1 <- simplify(net1, remove.loops = T)
        net1 <- as.undirected(net1)
        plot(net1,layout=layout.fruchterman.reingold(net1),vertex.label=NA,vertex.color=test$color,vertex.size=test$percent*10000)
        title(main=paste0(files_shortA,"len:",j,"Aas"))
      } else {
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      }
    }
  }
  dev.off()
}






