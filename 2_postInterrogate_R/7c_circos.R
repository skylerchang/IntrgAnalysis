library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(circlize)


t<-read_rds('../../Data/Clntab_RDS/clntab_vAndJ.rds')
datalist<-t[[1]]
files_short<-t[[2]]

outpath<-'../../Results/Circos/NumericalOrder/'



#create table combining the data from all samples
s.all<-NULL
for (i in 1:length(datalist)){
  s.all<-bind_rows(s.all,datalist[[i]])
}
s.all<-list(s.all)
datalist<-c(s.all,datalist)

#add file name for 'all' to filename vector
files_short<-c('allSamples',files_short)

gap<-vector()

for (i in 1:length(datalist)){
  s<-datalist[[i]]
  
  #s<-s[sample(nrow(s),1000),]
  
  #eliminate lines with ambiguous V or J call (i.e. lines containing '=')
  s<-s[!grepl("=",s$vGene),]
  s<-s[!grepl("=",s$jGene),]
  #eliminate NA recods
  s<-s[rowSums(is.na(s))==0,]
  s.all<-bind_rows(s.all,s)

  #subset according to locus
  aa<-s[s$vChain=='TRA' & s$jChain=='TRA',]
  ad<-s[s$vChain=='TRA' & s$jChain=='TRD',]
  bb<-s[s$vChain=='TRB' & s$jChain=='TRB',]
  da<-s[s$vChain=='TRD' & s$jChain=='TRA',]
  dd<-s[s$vChain=='TRD' & s$jChain=='TRD',]
  gg<-s[s$vChain=='TRG' & s$jChain=='TRG',]
  cc<-s[(s$vChain=='TRA' | s$vChain=='TRD') & (s$jChain=='TRA' | s$jChain=='TRD'),]
  
  #select vGene, jGene 7 size
  aa<-aa[,c(1,2,5)]
  ad<-ad[,c(1,2,5)]
  bb<-bb[,c(1,2,5)]
  da<-da[,c(1,2,5)]
  dd<-dd[,c(1,2,5)]
  gg<-gg[,c(1,2,5)]
  cc<-cc[,c(1,2,5)]
  
  #limit number of lines for easier plotting
#  aa<-aa[sample(nrow(aa),500),]
#  bb<-bb[sample(nrow(bb),500),]
#  dd<-dd[sample(nrow(dd),500),]
#  gg<-gg[sample(nrow(gg),500),]
#  cc<-cc[sample(nrow(gg),500),]
  
  #sum up rows with same v-j combination
  library(plyr)
  aa<-ddply(aa,c("vGene","jGene"), numcolwise(sum))
  bb<-ddply(bb,c("vGene","jGene"), numcolwise(sum))
  dd<-ddply(dd,c("vGene","jGene"), numcolwise(sum))
  gg<-ddply(gg,c("vGene","jGene"), numcolwise(sum))
  cc<-ddply(cc,c("vGene","jGene"), numcolwise(sum))
  detach("package:plyr",unload=T)
  
  #eliminate locus abbreviations
#  bb$vGene<-sub("TR[ABDG]","",bb$vGene)
#  bb$jGene<-sub("TR[ABDG]","",bb$jGene)
#  gg$vGene<-sub("TR[ABDG]","",gg$vGene)
#  gg$jGene<-sub("TR[ABDG]","",gg$jGene)
  


  #================== TRB ================== 
  #define order
  order.b<-c('TRBV1','TRBV4-1','TRBV2','TRBV3','TRBV4-2','TRBV5-1','TRBV5-2','TRBV6','TRBV7-1','TRBV8','TRBV5-3','TRBV7-2','TRBV5-4','TRBV10','TRBV11','TRBV12-1','TRBV12-2','TRBV15','TRBV16','TRBV18','TRBV19','TRBV20','TRBV21','TRBV22','TRBV23','TRBV24','TRBV25','TRBV26','TRBV27','TRBV28','TRBV29','TRBJ1-1','TRBJ1-2','TRBJ1-3','TRBJ1-4','TRBJ1-5','TRBJ2-1','TRBJ2-2','TRBJ2-3','TRBJ2-4','TRBJ2-5','TRBV30')
  length(order.b)
  
  grid.colors.b<-c(rep('#DF6F32',31),rep('grey',5),rep('grey',5),'#DF6F32')
  length(grid.colors.b)
  
  #define gaps dynamically, i.e. based on the number of genes present
  allUniqueVsAndJs<-c(unique(bb$vGene),unique(bb$jGene))
  length(allUniqueVsAndJs)
  
  exist<-order.b %in% allUniqueVsAndJs
  vs1<-sum(exist[1:31]==T)       #first block of Vs
  js1<-sum(exist[32:36]==T)      #first block of Js
  js2<-sum(exist[37:41]==T)      #second block of Js
  vs2<-sum(exist[42:42]==T)      #last V (TRBV30 in inverted orientation)
  #default gap
  #g<-c(rep(1,(vs1-1)),10,rep(1,(js1-1)),10,rep(1,(js2-1)),rep(10,vs2),10)
  length(gap)
  
  gap<-c(rep(1,(vs1-3)),10,rep(1,(js1-1)),10,rep(1,(js2-1)),rep(10,vs2),10) #adjusted gap 1
  
  grid.colors.b<-grid.colors.b[exist]
  length(grid.colors.b)

  pdf(paste0(outpath,files_short[i],'_trb.pdf'))
  circos.clear()
  circos.par(start.degree = 90)
  circos.par(gap.after=gap)
  chordDiagram(bb,annotationTrack = "grid", 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
  #chordDiagram(bb,order=order.b,grid.col=grid.colors.b,annotationTrack = "grid", 
   #            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.4)
  }, bg.border = NA)
  dev.off()
  
  #================== TRG ================== 
  order.g<-c('TRGV2-1','TRGJ1-1','TRGJ1-2','TRGV2-2','TRGJ2-1','TRGJ2-2','TRGV2-3','TRGJ3-1','TRGJ3-2','TRGJ4-3','TRGJ4-2','TRGJ4-1','TRGJ5-1','TRGV5-3','TRGV7-1','TRGV6-1','TRGV2-4','TRGJ6-1','TRGJ6-2')
  
  v2<-'red'
  v5<-'yellow'
  v7<-'blue'
  v6<-'green'
  #'#DF6F32'
  #'#CFA42D'
  
  grid.v1<-'black'
  grid.j1<-'light grey'
  
  grid.colors.g<-c(grid.v1,grid.j1,grid.j1,grid.v1,grid.j1,grid.j1,grid.v1,grid.j1,grid.j1,grid.j1,grid.j1,grid.j1,grid.j1,grid.v1,grid.v1,grid.v1,grid.v1,grid.j1,grid.j1)
  
  link.colors.g<-c(v2,'grey','grey',v2,'grey','grey',v2,'grey','grey','grey','grey','grey','grey',v5,v7,v6,v2,'grey','grey')
  
  
  #define gaps dynamically, i.e. based on the number of genes present
  allUniqueVsAndJs<-c(unique(gg$vGene),unique(gg$jGene))
  length(order.g)
  length(allUniqueVsAndJs)
  
  exist<-order.g %in% allUniqueVsAndJs
  cassette1<-sum(exist[1:3]==T)
  cassette2<-sum(exist[4:6]==T)
  cassette3<-sum(exist[7:9]==T)
  cassette4<-sum(exist[10:12]==T)
  cassette5<-sum(exist[13:16]==T)
  cassette6<-sum(exist[17:19]==T)
  gap<-c(rep(1,(cassette1-1)),10,rep(1,(cassette2-1)),10,rep(1,(cassette3-1)),10,rep(1,(cassette4-1)),10,rep(1,(cassette5-1)),10,rep(1,(cassette6-1)),10)
  length(gap)
  length(allUniqueVsAndJs)
  grid.colors.g<-grid.colors.g[exist]
#  link.colors.g<-link.colors.g[exist]
  
  pdf(paste0(outpath,files_short[i],'_trg.pdf'))
  circos.clear()
  circos.par(start.degree = 90)
  #circos.par(gap.after = gap)
  #m<-chordDiagram(gg,annotationTrack = "grid",
   #            link.lwd = 1, link.lty = 1, link.border = "black",
    #           preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
  m<-chordDiagram(gg,order=order.g,grid.col=grid.colors.g,col=link.colors.g,annotationTrack = "grid",
                  link.lwd = 1, link.lty = 1, link.border = "black",
                  preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.4)
  }, bg.border = NA)
  dev.off()
  
  
  #================= TRA/TRD ===================
  order.c<-c('TRAV2','TRAV3','TRAV8-1','TRAV4','TRAV43-1','TRAV5','TRAV6','TRAV9-1','TRAV8-2','TRAV9-2','TRAV8-3','TRAV43-2','TRAV14-1','TRAV13-1','TRAV9-3','TRAV14-2','TRAV9-4','TRAV43-3','TRAV14-3','TRAV9-5','TRAV8-4','TRAV43-4','TRAV13-2','TRAV9-6',"TRAV43-5",'TRAV43-6','TRAV13-3','TRAV14-5','TRAV9-7','TRAV10','TRAV12','TRAV8-5','TRAV16','TRAV17','TRAV18','TRAV19','TRAV20','TRAV21','TRAV8-6','TRAV22','TRAV23','TRDV1','TRAV24','TRAV25','TRAV8-7','TRAV27','TRAV28','TRAV29','TRDV2','TRAV26','TRAV34','TRAV35','TRAV36','TRAV37','TRAV38-1','TRAV38-2','TRAV40','TRAV41','TRDV5-1','TRDV5-2','TRDV5-3','TRDV5-4','TRDV5-5','TRDV5-6','TRDV6','TRDV4','TRDV3','TRDJ1','TRDJ4','TRDJ3','TRDJ2','TRAJ61','TRAJ65','TRAJ60','TRAJ59','TRAJ58','TRAJ57','TRAJ56','TRAJ54','TRAJ53','TRAJ64','TRAJ52','TRAJ51','TRAJ50','TRAJ49','TRAJ48','TRAJ47','TRAJ45','TRAJ44','TRAJ43','TRAJ42','TRAJ41','TRAJ40','TRAJ39','TRAJ38','TRAJ37','TRAJ36','TRAJ35','TRAJ34','TRAJ33','TRAJ32','TRAJ31','TRAJ30','TRAJ29','TRAJ28','TRAJ27','TRAJ26','TRAJ25','TRAJ24','TRAJ23','TRAJ22','TRAJ21','TRAJ20','TRAJ19','TRAJ18','TRAJ17','TRAJ16','TRAJ15','TRAJ63','TRAJ14','TRAJ13','TRAJ12','TRAJ11','TRAJ10','TRAJ9','TRAJ8','TRAJ7','TRAJ6','TRAJ5','TRAJ62','TRAJ4','TRAJ3','TRAJ2','TRAJ1')
  
  length(order.c)
  
  
  fc <- colorRampPalette(c("lightgreen", "darkgreen"))
  plot(rep(1, 10),col = fc(48), pch = 19, cex = 3)
  
  
  grid.colors.c<-c(fc(48),        #TRAVs 1-47
                   'gold',                          #TRDV 1
                   rep('cornflowerblue',9),         #TRAV 48-56
                   rep('gold',9),                   #TRDV 2-10
                   rep('lightgrey',4),              #TRDJ 1-4
                   rep('darkgrey',63))              #TRAJ 1-63
  length(grid.colors.c)
  
  #define gaps dynamically, i.e. based on the number of genes present
  allUniqueVsAndJs<-c(unique(cc$vGene),unique(cc$jGene))
  length(order.c)
  length(allUniqueVsAndJs)
  
  exist<-order.c %in% allUniqueVsAndJs
  v<-sum(exist[1:67]==T)
  trdj<-sum(exist[68:72]==T)
  traj<-sum(exist[73:134]==T)
  gap<-c(rep(1,(v-1)),10,rep(1,(trdj-1)),10,rep(1,(traj-1)),10)
  
  gap<-c(rep(1,(v-3)),10,rep(1,(trdj-1)),10,rep(1,(traj-2)),10) #gap 'all'
  
  length(gap)
  length(allUniqueVsAndJs)
  grid.colors.c<-grid.colors.c[exist]
  length(grid.colors.c)
  
  pdf(paste0(outpath,files_short[i],'_tra-trd.pdf'))
  circos.clear()
  circos.par(start.degree = 90)
  circos.par(gap.after = gap)
  chordDiagram(cc,annotationTrack = "grid", 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
  chordDiagram(cc,order=order.c,grid.col=grid.colors.c,annotationTrack = "grid", 
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.4)
  }, bg.border = NA)
  dev.off()
}


#===== test example ======
s<-tibble(vGene=c('TRGV1','TRBV1','TRAV1','TRDV1','TRAV1','TRDV1'),
          jGene=c('TRGJ1','TRBJ1','TRAJ1','TRDJ1','TRDJ1','TRAJ1'),
          aaSeq=rep('asdf',6),
          aaLength=rep('asdf',6),
          size=c(1,1,1,1,1,1),
          vChain=c('TRG','TRB','TRA','TRD','TRA','TRD'),
          jChain=c('TRG','TRB','TRA','TRD','TRD','TRA'))

order.c<-c('TRAV1','TRDV1','TRDJ1','TRAJ1')

#set grid colors
grid.colors<-c(rep('grey',length(cc$vGene)),rep('black',length(cc$jGene)))

#sum up rows with same v-j combination
library(plyr)
cc<-ddply(cc,c("vGene","jGene"), numcolwise(sum))
detach("package:plyr",unload=T)

circos.clear()
circos.par(start.degree = 90)
chordDiagram(cc,order=order.c)
