setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_2018-03-14_v5/')

#*** create table tMain **********
#import read & uniread number tables
reads<-read.table('Results/clntab_ALL-reads.txt', header=T)
unireads<-read.table('Results/clntab_ALL-unireads.txt', header=T)
#combine the 1st three columns of each table into a new table
main<-cbind(reads[,1:3],unireads[,1:3])
#change the column names to indicate read & uniread category
colnames(main)<-c("TRB.read","TRD.read","TRG.read","TRB.uniread","TRD.uniread","TRG.uniread")

main$TRB.div<-main$TRB.read/main$TRB.uniread
main$TRD.div<-main$TRD.read/main$TRD.uniread
main$TRG.div<-main$TRG.read/main$TRG.uniread

#import sample information; I saved the xlsx file as tab deliminted file
details<-read.table('Results/sampleDetails.txt', header=T)
main<-cbind(main,details)


library(rgl)
organ= as.factor(main$organ)  #Levels: Li Ln SI Sk Sp
colbyorgan = c("purple3", "green4", "salmon2", "blue3", "red3")[organ] # Li → purple, Ln → green, SI → salmon, Sk → blue3, Sp → red

open3d()
ReadPlot <- plot3d(main$TRB.read, main$TRD.read, main$TRG.read, type = "p", size = 8, col = colbyorgan, xlab="TRB", ylab="TRD", zlab="TRG", main="Reads") #plot graph1

open3d() # to open multiple rgl windows
UniReadPlot<-plot3d(main$TRB.uniread, main$TRD.uniread, main$TRG.uniread, type = "p", size = 8, col = colbyorgan, xlab="TRB", ylab="TRD", zlab="TRG", main="UniReads") #plot graph 2

open3d()
DivPlot <- plot3d(main$TRB.div, main$TRD.div, main$TRG.div, type = "p", size = 6, col = colbyorgan, xlab="TRBdiv", ylab="TRDdiv", zlab="TRGdiv", main="Reads/Unireads") # plot diversity graph

test<-main[,1:9]
test.pca<-prcomp(test, center = TRUE, scale. = TRUE) 
plot(test.pca)
summary(test.pca)
library(ggbiplot)
g <- ggbiplot(test.pca, obs.scale = 1, var.scale = 1, groups = main$organ, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
