library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(seqinr)
library(plyr)
library(igraph)
library(stringdist)

t<-read_rds('RDS/clntab_vAndJ.rds')
outpath<-'OUT/cdrAaComparison'
loci<-c("TRB","IGH")
#Top n clones that are being displayed in separate colors (all other clones are grey)
aanum<-100
#==========================================================
#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]
print(files_short)
my.sample.1<- readline(prompt = "Enter sample1 name: ")
my.sample.2<- readline(prompt = "Enter sample2 name: ")
sampleID.1<-which(grepl(my.sample.1,t[[2]]))
sampleID.2<-which(grepl(my.sample.2,t[[2]]))
compare_list<-list(datalist[[sampleID.1]],datalist[[sampleID.2]])
name_list<-list(my.sample.1,my.sample.2)

for (i in 1:length(loci)) {
      pdf(paste0(outFolder,"/cln2igr_",loci[i],"_loci.pdf"),width=9,height=11)
      par(mfrow=c(4,4),mar=c(1,1,1,1))
      for (j in seq(from=8,to=23,by=1))  {
        for (k in 1:length(compare_list)) {
          d<-compare_list[[k]][,c("vGene","jGene","aaSeq","aaLength","size","vAndJchainSimplified")]
          dd<-d[d$vAndJchainSimplified==loci[i],c("aaSeq","size","aaLength")]
          dd<-dd[!is.na(dd$aaSeq),]
          ### need  to change ######
          files_shortA <- gsub(".*[AB]_","",name_list[[k]])
          tt2<-as_tibble(ddply(dd,.(aaSeq,aaLength),numcolwise(sum)))
          tt2<-tt2[order(-tt2$size),]
          tt2$aaLength[tt2$aaSeq=="noju"]<-"4"
          tt2<-tt2[complete.cases(tt2),]
          tt2$color<-"grey"
          tt2$color[1:100]<-color100
          #test$color[test$percent < 0.1] <-"grey"
          #test<-tt2[tt2$aaLength==j,c("aaSeq","size","color")]
          test<-tt2[tt2$aaLength=="13",c("aaSeq","size","color")]
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
