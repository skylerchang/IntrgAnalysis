library(RColorBrewer)
setwd('/Users/SKeller/Documents/Sequencing/Runs/k9MultiLoci3-73454401/Clntab_2018-04-25/Results')

for (type in c('reads','unireads')){
  #t<-read.table('clntab_ALL-unireads.txt', sep = " ", header=T, check.names = FALSE)
  t<-read.table(paste0('clntab_ALL-',type,'.txt'), sep = " ", header=T, check.names = FALSE)
  
  #subset samples by project
  t<-t[c(5,6,1,2,3,4),]   #Wylde Blue
#  t<-t[c(7:24),]          #Ultra II kit vs HB
  
  #shorten file name
  t$id<-sub(".*Data/.*/", "",t$sample)
  t$id<-sub("_.*", "",t$id)
  
  #subset IGH,TRB,TRD,TRG
  l<-t[,c(1,3,4,5)]
  
  #pdf("IGH-TRB-TRD-TRG-ratio_unireads.pdf")
  #png("IGH-TRB-TRD-TRG-ratio_unireads.png")
  pdf(paste0("IGH-TRB-TRD-TRG-ratio_",type,".pdf"))
#  png(paste0("IGH-TRB-TRD-TRG-ratio_",type,".png"))
  barplot(t(as.matrix(l)),names.arg=t$id,las=2,col=brewer.pal(ncol(l),"Set3"),main=paste0(type,' per locus'))
  legend("topleft",legend=colnames(l),fill=brewer.pal(ncol(l),"Set3"),cex=0.8)
  dev.off()
}

