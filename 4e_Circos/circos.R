setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_test/')

library(circlize)

dirs <- list.dirs(path = 'Data/', full.names=T, recursive=F)
for (i in 1:length(dirs)) {
  print(dirs[i])
  
  datalist <- list()
  
  files <- list.files(path = dirs[i], pattern="TR[B]_v-j-size", full.names=T, recursive=F)
  for (j in 1:length(files)) {
    print(files[j])
    t <- read.table(files[j], sep="\t", colClasses = c(rep("character",2),"integer"),header = TRUE,fill = TRUE,comment.char="")
    t<- head(t[order(-t$size),],100)
    
    circos.clear()
    chordDiagram(t)
  }
}

