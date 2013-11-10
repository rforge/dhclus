withinSum<-function(Data,labels,Distance=FALSE){
Data<-t(Data)
if(!Distance)DMAT<-.Call("icav_sq_euc_metric",Data, c(dim(Data)[1],dim(Data)[2]),DUP=F,PACKAGE="icav")
else DMAT<-Data
if (ncol(DMAT)!=length(labels)) stop("Data length and labels length do not match")
W<-.Call("icav_gap_withinSum",DMAT,as.integer(labels),unique(labels),length(labels),DUP=F,PACKAGE="icav")
return(W)
}
