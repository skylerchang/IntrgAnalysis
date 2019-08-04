
library(plyr)
library(tidyverse)
library(gridExtra)
library(here)
library(openxlsx)
library(RColorBrewer)

setwd(here())

source("2_postInterrogate_R/functions.R")

#========== adjust the following variables ================
sampleNoCodesForFraction<-F        #T or F
run<-25
loci<-c('IGH','TRB','TRG')

t<-read_rds('../Data/clntab_RDS/clntab_vAndJ_filtered.rds')
outpath<-'../Results/GeneUsage/'

#==========================================================
dir.create(outpath)

#datalist contains one tibble per sample
datalist<-t[[1]]
#sample names are stored separately in a vector
files_short<-t[[2]]

data<-list()
factorLevels<-list()

#initialize data list that holds summary stats for groups
groups<-c('vSubgroup','vGene','jSubgroup','jGene')
for (j in 1:length(groups)){
  data[[j]]<-tibble()
}

for (i in 1:length(datalist)){
  print(files_short[[i]])
  
  #parse filtered clntab data 
  t<-datalist[[i]][,c("vGene","jGene","readCount","locus")]
  if (is.null(datalist[[i]])){
    print("no data")
    next
  }
  
  #determine subgroup
  findSubgroups<-function(x){
    x %>% 
      gsub("-[0-9]+","",.) %>%
      strsplit(split="=") %>%
      lapply(unique) %>%
      sapply(paste,collapse="=")
  }
  t$vSubgroup<-findSubgroups(t$vGene)
  t$jSubgroup<-findSubgroups(t$jGene)

  #group & summarize - function
  group<-function(x,gene){
    x %>%
      rename(group=!!gene) %>%
      group_by(locus) %>%
      mutate(readCountPerLocus=sum(readCount)) %>%
      mutate(percentage=readCount/readCountPerLocus*100) %>%
      group_by(locus,group) %>%
      summarise(percentage=sum(percentage),readCount=sum(readCount)) %>%
      mutate(filename=files_short[i])
  }
  #group, summarise & add to table that holds data from all files
  for (j in 1:length(groups)){
    t.grouped<-group(t,groups[j])
    data[[j]]<-bind_rows(data[[j]],t.grouped)
  }
}

#plot
print("plotting...")
for (i in 1:length(groups)){
  print(groups[i])
  d<-data[[i]]
  #determine overall rank based on frequency -> determines factor level for plotting
  factorLevels[[i]]<-d %>%
    group_by(locus,group) %>% 
    summarise(readCount=sum(readCount)) %>%
    arrange(locus,desc(readCount)) %>%
    pull(group)
  #modify factor level for 'group'  
  d$group<-factor(d$group,levels=factorLevels[[i]]) 
  
  #split filename
  d<-splitFilename(d,sampleNoCodesForFraction,run)
  #split data by locus to enable plots with individual scales
  s<-split(d,d$locus)
  
  #facet by submission
  types<-c('percentage','readCount')
  for (k in 1:length(types)){
    p<-list()
    for (l in 1:length(s)){
      p[[l]]<-ggplot(s[[l]],aes_string(types[k],"group",fill="submission"))+
        geom_point(pch=21)+
        facet_grid(submission~locus)+
        theme(legend.position="none")
    }
    pdf(paste0(outpath,"geneUsage_",groups[i],"_",types[k],".pdf"))
    grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)
    dev.off()
  }
}

