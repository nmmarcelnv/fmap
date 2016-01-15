args=commandArgs(trailingOnly = TRUE)
nb=as.integer(args[2])
x<-readBin(args[1],what="double",n=nb*nb)
xm<-as.matrix(x,c(nb,nb))

