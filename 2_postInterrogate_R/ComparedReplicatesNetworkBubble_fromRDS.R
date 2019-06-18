library(tidyverse)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(seqinr)
library(plyr)
library(igraph)
library(stringdist)
library(kader)
library(here)

#####assign 100 colors to bubbles#######
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


######Create Network plots for Aa bubbles########

setwd(here())
getwd()

t<-read_rds('../Data/Clntab_RDS/clntab_vAndJ.rds')
outpath<-'../Results/Networkplots/'
loci<-c("TRB","IGH")
#Top n clones that are being displayed in separate colors (all other clones are grey)
aanum<-100
#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

dir.create(outpath<-'../Results/Networkplots/',recursive=T)

###### Compare every 2 replicates of same sample#########
for (m in seq(from=1,to=length(datalist),by=2)) {
my.sample.1<- files_short[m]
my.sample.2<- files_short[m+1]
name_list<-list(my.sample.1,my.sample.2)
files_shortA<- gsub("-[AB].*","",my.sample.1)
full_names<- gsub("_S.","",my.sample.1)
compare_list<-list(datalist[[m]],datalist[[(m+1)]])
####### test different locus ########
  for (i in 1:length(loci)) {
      pdf((paste0(outpath,files_shortA,"-",loci[i],"_loci.pdf")),width=10,height=13)
      par(mfrow=c(4,4),mar=c(1,1,1,1),oma=c(0,0,2,0))
      ############ test on Amino acid length from 8 to 23########
      for (j in seq(from=8,to=23,by=1))  {
        for (k in 1:length(compare_list)) {
          d<-compare_list[[k]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]
          dd<-d[d$vAndJchainSimplified==loci[i],c("aaSeq","size","aaLength")]
          dd<-dd[!is.na(dd$aaSeq),]
          tt2<-as_tibble(ddply(dd,.(aaSeq,aaLength),numcolwise(sum)))
          tt2<-tt2[order(-tt2$size),]
          tt2$aaSeq[tt2$aaSeq=="noju"]<-"NA"
          tt2<-tt2[complete.cases(tt2),]
          if ( nrow(tt2) > 100 ) { 
          tt2$color<-"grey"
          tt2$color[1:100]<-col
         } else if (nrow(tt2) > 0 && nrow(tt2)<100 ) {
          tt2$color<-"grey"
          tt2$color[1:nrow(tt2)]<-col[1:nrow(tt2)]
        } else {
          break
        }  
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
            plot(net2,layout=layout.fruchterman.reingold(net2),vertex.label=NA,vertex.color=test2$color,
                 vertex.size=kader:::cuberoot(test2$percent*10e5),main=paste0("Amino acid length:",j,"Aas"))
          } else if (nrow(test) > 0 && nrow(test) < aanum) { 
            uniqueids <- as.character(test$aaSeq)
            hdids <- stringdistmatrix(uniqueids,uniqueids,method = "hamming")  
            rownames(hdids) <- uniqueids
            colnames(hdids) <- uniqueids
            links <- hdids
            net1 <- graph_from_adjacency_matrix(links<3)
            net1 <- simplify(net1, remove.loops = T)
            net1 <- as.undirected(net1)
            plot(net1,layout=layout.fruchterman.reingold(net1),vertex.label=NA,vertex.color=test$color,
                 vertex.size=kader:::cuberoot(test$percent*10e5), main=paste0("Amino acid length:",j,"Aas"))
          } else {
            plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
          }
        }
      }
      mtext((paste0("Network Bubble plots of ",full_names)), outer=TRUE, cex=1,side=3)
      dev.off()
      } 
       if (m > length(datalist)) {
      break
    }   
}
       
