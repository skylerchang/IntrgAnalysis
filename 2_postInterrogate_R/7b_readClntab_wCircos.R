library(tidyverse)
library(here)
library(seqinr)
library(RColorBrewer)
library(gridExtra)
library(circlize)

targetFolder<-"../../Data/Clntab/"
outFolder<-"../../Data/CountPlots/"

datalist<-list()

#list all files in Clntab folder, which have long suffix
files<-list.files(targetFolder,pattern ="*.clntab",recursive = T)
#truncate suffix
files_short<-basename(sub("_L001_R1_001.fastq.processed.junctioned.profiled.clntab","",files))


plotlist<-list()
datalist<-list()

#for (i in 1:length(files)){
for (i in 1:2){
  print(files[i])
  t<-read_tsv(paste0(targetFolder,files[i]))
  d<-subset(t, select=c("sequence.5-GENE","sequence.3-GENE","sequence.JUNCTION.aa seq","sequence.JUNCTION.aa seq.len","sequence.size"))
  colnames(d)<-c("vGene","jGene","aaSeq","aaLength","size")
  
  d$vChain<-NA
  d$vChain[grepl("TRA",d$vGene)]<-'TRA'
  d$vChain[grepl("TRB",d$vGene)]<-'TRB'
  d$vChain[grepl("TRD",d$vGene)]<-'TRD'
  d$vChain[grepl("TRG",d$vGene)]<-'TRG'
  
  d$jChain<-NA
  d$jChain[grepl("TRA",d$jGene)]<-'TRA'
  d$jChain[grepl("TRB",d$jGene)]<-'TRB'
  d$jChain[grepl("TRD",d$jGene)]<-'TRD'
  d$jChain[grepl("TRG",d$jGene)]<-'TRG'
  
  d$junction<-NA
  d$junction[d$vGene=='no5' & d$jGene=='no3']<-'neitherVnorJ'
  d$junction[d$vGene=='no5' & d$jGene!='no3']<-'jOnly'
  d$junction[d$vGene!='no5' & d$jGene=='no3']<-'vOnly'
  d$junction[d$vGene!='no5' & d$jGene!='no3']<-'vAndJ'
  
  d$vAndJchain<-paste(d$vChain,d$jChain,sep='-')
  
  d$vAndJchainSimplified<-NA
  d$vAndJchainSimplified[d$vAndJchain=='NA-NA' | d$vAndJchain=='TRD-TRA']<-'Other'
  d$vAndJchainSimplified[grepl('TRA',d$vAndJchain)]<-'TRA'
  d$vAndJchainSimplified[grepl('TRB',d$vAndJchain)]<-'TRB'
  d$vAndJchainSimplified[grepl('TRD',d$vAndJchain)]<-'TRD'
  d$vAndJchainSimplified[grepl('TRG',d$vAndJchain)]<-'TRG'

  d$junction<-factor(d$junction,levels=c('vAndJ','jOnly','neitherVnorJ'))
  
  #reduce number of lines to speed up plotting
  dd<-d[sample(nrow(d),1000),]
  
  p<-ggplot(dd,aes(junction,size,fill=vAndJchainSimplified))+geom_col()+scale_fill_brewer(palette='Set3')
  
  plotlist[[i]]<-ggplotGrob(p)

  datalist[[i]]<-d[d$junction=='vAndJ',]
}

marrangeGrob(plotlist,nrow=2,ncol=2)

s<-datalist[[1]]

#s<-s[sample(nrow(s),1000),]

#eliminate lines with ambiguous V or J call (i.e. lines containing '=')
s<-s[!grepl("=",s$vGene),]
s<-s[!grepl("=",s$jGene),]
#eliminate NA recods
s<-s[rowSums(is.na(s))==0,]
#select cols vGene, jGene, size

#subset according to locus
aa<-s[s$vChain=='TRA' & s$jChain=='TRA',]
ad<-s[s$vChain=='TRA' & s$jChain=='TRD',]
bb<-s[s$vChain=='TRB' & s$jChain=='TRB',]
da<-s[s$vChain=='TRD' & s$jChain=='TRA',]
dd<-s[s$vChain=='TRD' & s$jChain=='TRD',]
gg<-s[s$vChain=='TRG' & s$jChain=='TRG',]
cc<-s[(s$vChain=='TRA' | s$vChain=='TRD') & (s$jChain=='TRA' | s$jChain=='TRD'),]

aa<-aa[,c(1,2,5)]
ad<-ad[,c(1,2,5)]
bb<-bb[,c(1,2,5)]
da<-da[,c(1,2,5)]
dd<-dd[,c(1,2,5)]
gg<-gg[,c(1,2,5)]



#limit number of lines for easier plotting
aa<-aa[sample(nrow(aa),500),]
bb<-bb[sample(nrow(bb),500),]
dd<-dd[sample(nrow(dd),500),]
gg<-gg[sample(nrow(gg),500),]
cc<-cc[sample(nrow(gg),500),]

#sum up rows with same v-j combination
library(plyr)
aa<-ddply(aa,c("vGene","jGene"), numcolwise(sum))
bb<-ddply(bb,c("vGene","jGene"), numcolwise(sum))
dd<-ddply(dd,c("vGene","jGene"), numcolwise(sum))
gg<-ddply(gg,c("vGene","jGene"), numcolwise(sum))
cc<-ddply(cc,c("vGene","jGene"), numcolwise(sum))
detach("package:plyr",unload=T)

#eliminate locus abbreviations
bb$vGene<-sub("TR[ABDG]","",bb$vGene)
bb$jGene<-sub("TR[ABDG]","",bb$jGene)
gg$vGene<-sub("TR[ABDG]","",gg$vGene)
gg$jGene<-sub("TR[ABDG]","",gg$jGene)

#================== TRB ================== 
#define order
order.b<-c('TRBV1','TRBV4-1','TRBV2','TRBV3','TRBV4-2','TRBV5-2','TRBV6','TRBV7-1','TRBV8','TRBV5-3','TRBV7-2','TRBV5-4','TRBV10','TRBV11','TRBV12-1','TRBV12-2','TRBV15','TRBV16','TRBV18','TRBV19','TRBV20','TRBV21','TRBV22','TRBV23','TRBV24','TRBV25','TRBV26','TRBV27','TRBV28','TRBV29','TRBJ1-1','TRBJ1-2','TRBJ1-3','TRBJ1-4','TRBJ1-5','TRBJ2-1','TRBJ2-2','TRBJ2-3','TRBJ2-4','TRBJ2-5','TRBV30')

grid.colors.b<-c(rep('#DF6F32',30),rep('grey',5),rep('grey',5),'#DF6F32')

#reduce number of lines for faster plotting
#define gaps dynamically, i.e. based on the number of genes present
allUniqueVsAndJs<-c(unique(bb$vGene),unique(bb$jGene))
length(order.b)
length(grid.colors.b)
length(allUniqueVsAndJs)

exist<-order.b %in% allUniqueVsAndJs
vs1<-sum(exist[1:30]==T)       #first block of Vs
js1<-sum(exist[31:35]==T)      #first block of Js
js2<-sum(exist[36:40]==T)      #second block of Js
vs2<-sum(exist[41:41]==T)      #last V (TRBV30 in inverted orientation)
gap<-c(rep(1,(vs1-1)),10,rep(1,(js1-1)),10,rep(1,(js2-1)),10,rep(10,vs2))
length(gap)
length(allUniqueVsAndJs)
grid.colors.b<-grid.colors.b[exist]

pdf('feTRB.pdf')
circos.clear()
circos.par(start.degree = 90)
circos.par(gap.after=gap)
chordDiagram(bb,order=order.b,grid.col=grid.colors.b,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.4)
}, bg.border = NA)
dev.off()

#================== TRG ================== 
order.g<-c('TRGV2-1','TRGJ1-1','TRGJ1-2','TRGV2-2','TRGJ2-1','TRGJ2-2','TRGV2-3','TRGJ3-1','TRGJ3-2','TRGJ4-3','TRGJ4-2','TRGJ4-1','TRGJ5-1','TRGV5-3','TRGV7-1','TRGV6-1','TRGV2-4','TRGJ6-1','TRGJ6-2')

v2<-'#E41A1C'
v5<-'#3C899E'
v7<-'#49A75A'
v6<-'#D28AB0'
#'#DF6F32'
#'#CFA42D'

grid.colors.g<-c(v2,'grey','grey',v2,'grey','grey',v2,'grey','grey','grey','grey','grey','grey',v5,v7,v6,v2,'grey','grey')


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

pdf('feTRG.pdf')
circos.clear()
circos.par(start.degree = 90)
circos.par(gap.after = gap)
chordDiagram(gg,order=order.g,grid.col=grid.colors.g,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.4)
}, bg.border = NA)
dev.off()


#================= TRA/TRD ===================
order.c<-c('TRAV2','TRAV3','TRAV8-1','TRAV4','TRAV43-1','TRAV5','TRAV6','TRAV9-1','TRAV8-2','TRAV9-2','TRAV8-3','TRAV43-2','TRAV14-1','TRAV13-1','TRAV9-3','TRAV14-2','TRAV9-4','TRAV43-3','TRAV14-3','TRAV9-5','TRAV8-4','TRAV43-4','TRAV13-2','TRAV9-6','TRAV43-6','TRAV13-3','TRAV14-5','TRAV9-7','TRAV10','TRAV12','TRAV8-5','TRAV16','TRAV17','TRAV18','TRAV19','TRAV20','TRAV21','TRAV8-6','TRAV22','TRAV23','TRAV24','TRAV25','TRAV8-7','TRAV27','TRAV28','TRAV29','TRDV2','TRAV26','TRAV34','TRAV35','TRAV36','TRAV37','TRAV38-1','TRAV38-2','TRAV40','TRAV41','TRDV5-1','TRDV5-2','TRDV5-3','TRDV5-4','TRDV5-5','TRDV5-6','TRDV6','TRDV4','TRDV3','TRDJ1','TRDJ4','TRDJ3','TRDJ2','TRAJ61','TRAJ65','TRAJ60','TRAJ59','TRAJ58','TRAJ57','TRAJ56','TRAJ54','TRAJ53','TRAJ64','TRAJ52','TRAJ51','TRAJ50','TRAJ49','TRAJ48','TRAJ47','TRAJ45','TRAJ44','TRAJ43','TRAJ42','TRAJ41','TRAJ40','TRAJ39','TRAJ38','TRAJ37','TRAJ36','TRAJ35','TRAJ34','TRAJ33','TRAJ32','TRAJ31','TRAJ30','TRAJ29','TRAJ28','TRAJ27','TRAJ26','TRAJ25','TRAJ24','TRAJ23','TRAJ22','TRAJ21','TRAJ20','TRAJ19','TRAJ18','TRAJ17','TRAJ16','TRAJ15','TRAJ63','TRAJ14','TRAJ13','TRAJ12','TRAJ11','TRAJ10','TRAJ9','TRAJ8','TRAJ7','TRAJ6','TRAJ5','TRAJ62','TRAJ4','TRAJ3','TRAJ2','TRAJ1')

length(order.c)

grid.colors.c<-c(rep('cornflowerblue',46),        #TRAVs 1-46
                 'gold',                          #TRDV 1
                 rep('cornflowerblue',9),         #TRAV 47-58
                 rep('gold',9),                   #TRDV 2-10
                 rep('lightgrey',4),              #TRDJ 1-4
                 rep('darkgrey',63))              #TRAJ 1-63
length(grid.colors.c)

#define gaps dynamically, i.e. based on the number of genes present
allUniqueVsAndJs<-c(unique(cc$vGene),unique(cc$jGene))
length(order.c)
length(allUniqueVsAndJs)

exist<-order.c %in% allUniqueVsAndJs
v<-sum(exist[1:65]==T)
trdj<-sum(exist[66:70]==T)
traj<-sum(exist[71:132]==T)
gap<-c(rep(1,(v-1)),10,rep(1,(trdj-1)),10,rep(1,(traj-1)),10)
length(gap)
length(allUniqueVsAndJs)
#grid.colors.c<-grid.colors.c[exist]


circos.clear()
circos.par(start.degree = 90)
circos.par(gap.after = gap)

#chordDiagram(cc[,1:3],order=order.c,grid.col=grid.colors.c)

pdf("feTRA-TRD.pdf")
chordDiagram(cc[,1:3],order=order.c,grid.col=grid.colors.c,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cc))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.4)
}, bg.border = NA)
dev.off()

title("feTRA/TRD",cex = 0.8)

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


#======================================================





#set link colors
colorCount<-length(unique(s$vGene))+length(unique(s$jGene))
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))
fill=getPalette(colorCount)

#set link colors - ramp
col_fun = colorRamp2(range(t$size), c("#FFEEEE", "#FF0000"), transparency = 0.5)



for (i in 1:length(loci)){
  c<-loci[i]
  #draw plot
  circos.clear()
  #circos.par(gap.after=gaps)
  #circos.text=(facing = "clockwise")
  chordDiagram(c)
}


c


seqs<-as_tibble(seqs)
require(plyr)
seqs$concat<-paste(seqs$file,seqs$aaSeq)
summary<-ddply(seqs,.(aaSeq),summarise,freq=length(aaSeq))
summary<-head(summary[order(-summary$freq),],n=10)
topTen<-summary$aaSeq
seqs$aaSeq<-factor(seqs$aaSeq)
seqs$file<-factor(seqs$file)

seqsTop<-seqs[seqs$aaSeq %in% topTen,]


colorCount<-nlevels(seqs$file)
getPalette<-colorRampPalette(brewer.pal(9,'Set1'))
pdf(paste0(outFolder,"topTenClones_cloneVsAbundance.pdf"))
ggplot(seqsTop,aes(aaSeq,fill=file))+geom_bar()+coord_flip()+scale_fill_manual(values = getPalette(colorCount))+guides(fill=guide_legend(title = 'Samples',ncol = 1))
dev.off()

pdf(paste0(outFolder,"topTenClones_sampleVsClone.pdf"))
ggplot(seqsTop,aes(file,aaSeq))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



plotlistLoci<-list()
plotlistTRB<-list()
plotlistIGH<-list()

for (i in 1:length(datalist)){
  t<-datalist[[i]]
  p<-ggplot(t,aes(chain, fill=(validCdr)))+geom_bar()+ggtitle(files_short[i])+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('A','B','D','G','H','K','NA'))
  plotlistLoci[[i]]<-ggplotGrob(p)
  
  t.b<-t[t$validCdr==T & t$chain=='B',]
  p<-plotlistTRB[[i]]<-ggplot(t.b,aes(aaLength,fill=(aaSeq)))+geom_bar()+guides(fill=F)+xlim(4,38)+ggtitle(files_short[i])+theme(plot.title = element_text(hjust = 0.5,size=10))
  plotlistTRB[[i]]<-ggplotGrob(p)
  
  t.h<-t[t$validCdr==T & t$chain=='H',]
  p<-plotlistIGH[[i]]<-ggplot(t.h,aes(aaLength,fill=(aaSeq)))+geom_bar()+guides(fill=F)+xlim(4,38)+ggtitle(files_short[i])+theme(plot.title = element_text(hjust = 0.5,size=10))
  plotlistIGH[[i]]<-ggplotGrob(p)
}


require(gridExtra)
pdf(paste0(outFolder,'clntab_byLocus.pdf'))
marrangeGrob(plotlistLoci, nrow=2, ncol=2)
dev.off()

pdf(paste0(outFolder,'clntab_trbByAaLength.pdf'))
print(marrangeGrob(plotlistTRB, nrow=4, ncol=4))
dev.off()

pdf(paste0(outFolder,'clntab_ighByAaLength.pdf'))
print(marrangeGrob(plotlistIGH, nrow=4, ncol=4))
dev.off()

ggplot(seqs,aes(aaSeq))+geom_bar()


#extract nonsense reads
i<-1
t<-read_tsv(paste0(targetFolder,files[i]))
tt<-subset(t, select=c("sequence.5-GENE","sequence.3-GENE","sequence.JUNCTION.aa seq","sequence.JUNCTION.aa seq.len","sequence.chain","sequence.nt seq"))
str(tt)  
colnames(tt)<-c("v","j","aaSeq","aaLength","chain","nt")
trbj16<-tt$nt[tt$j=='TRBJ1-6']
noju<-tt$nt[tt$aaSeq=='noju']
write.fasta(as.list(noju[1:100]),rep("trbj16",100),'trbj16_100.fa')
