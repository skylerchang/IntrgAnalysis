#functions that are used by multiple scripts

splitFilename<-function(x,y){
  splitFilename<-strsplit(x$filename,"_")
  x$id<-laply(splitFilename, '[[', 1)
  x$ownerPatient<-laply(splitFilename, '[[', 2)
  x$idNumber<-laply(splitFilename, '[[', 3)
  x$pcr<-sub("C[0-9]+P[0-9]+C[0-9]+","",x$id)
  x$pcr<-sub("C[0-9]+L[0-9]+P[0-9]+","",x$pcr)
  x$sample<-sub("-*D[0-9]+P[0-9]+$","",x$pcr)
  x$submission<-sub("-[0-9a-zA-Z]+$","",x$sample)
  x$replicate<-ifelse(grepl("D[0-9]+P[0-9]*[13579]$",x$pcr),"rep1","rep2")
  if(y==T){
    x$fraction<-ifelse(grepl("-1$",x$sample),"PBMCs","Plasma")
  }else{
    x$fraction<-rep("fraction1",length(sample))
  }
  x
}