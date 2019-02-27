
# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)

#*********** test2 example **********
setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_test2/Data/')
vertices<-read.table('vertices.txt',header=T,colClasses=c('character','numeric'))
edges<-read.table('edges.txt',header=T,colClasses=c('character','character'))
# ***************************

#*********** test example **********
setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_test/Data/')
t <- read.table("1A_S23/test_top200.clntab",header=F,fill = TRUE,colClasses=c('character','character','character','numeric'))
colnames(t)<-c("aaSeq","length","functionality","size")

#***** edges ******
edges<-t[,c(2,1)]
colnames(edges)<-c('from','to')
u<-sort(unique(edges$from))
o<-rep('origin',length(u))
first<-data.frame('from'=o,'to'=u)
edges<-rbind(edges,first)

#***** vertices ******
vertices<-t[,c(1,4)]
colnames(vertices)<-c('name','size')
s<-rep(1,length(u))
firstV<-data.frame('name'=u,'size'=s)
vertices<-rbind(vertices,firstV)
origin<-data.frame('name'='origin','size'=1)
vertices<-rbind(vertices,origin)

#******* plot *********
mygraph <- graph_from_data_frame( edges, vertices=vertices)

p=ggraph(mygraph, layout = 'circlepack', weight="size") + 
  geom_node_circle(aes(fill = depth)) +
  theme_void() + 
  theme(legend.position="FALSE")
p

# Adjust color palette:
#p + scale_fill_viridis()
#p + scale_fill_distiller(palette = "RdPu") 


#> head(flare$edges)
#from                                           to
#1 flare.analytics.cluster flare.analytics.cluster.AgglomerativeCluster

#> head(flare$vertices)
#name size             shortName
#1 flare.analytics.cluster.AgglomerativeCluster 3938  AgglomerativeCluster

