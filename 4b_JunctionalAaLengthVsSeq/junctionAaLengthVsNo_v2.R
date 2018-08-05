setwd('/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_test/Data/')
library(ggplot2)

t <- read.table("1A_S23/1A_S23_TRB_col15-21-23-27.clntab",header=F,fill = TRUE)
colnames(t)<-c("aaSeq","length","functionality","size")

t$aaSeq <- reorder(t$aaSeq, t$size)
t<-head(t, n=100)

reds<-c(hsv(.01,.9,.9),hsv(.01,.9,.85),hsv(.01,.9,.8),hsv(.01,.9,.75),hsv(.01,.9,.7))
greens<-c(hsv(.4,.9,.9),hsv(.4,.9,.85),hsv(.4,.9,.8),hsv(.4,.9,.75),hsv(.4,.9,.7))
t$color <- ifelse(t$functionality == "productive", sample(greens), sample(reds))

ggplot(data = t, aes(x = length, y = size, fill = aaSeq)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=t$color, guide = guide_legend(reverse=TRUE, limits='CSVLLCQLPGGR#YERYF')) 



ggplot(data = t, aes(x = length, y = size, fill = aaSeq)) + 
  geom_bar(stat = "identity") + scale_fill_brewer(palette="Dark2", guide = guide_legend(reverse=TRUE))

ggplot(data = t, aes(x = length, y = size, fill = aaSeq)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values=t$color, guide = guide_legend(reverse=TRUE)) 
