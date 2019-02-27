setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_2018-03-14_v5/Results/')

t<-read.table('clntab_ALL-unireads.txt',header=T)

library(rgl)
plot3d(t$TRB,t$TRD,t$TRG,size=5,col='darkolivegreen3')