eucmetric<-function(data){

euc<-.Call("icav_euc_metric",data, c(dim(data)[1],dim(data)[2]),DUP=F)
#par[0] rows (dimension of Data), par[1] cols (amount of Data)
return(euc)
}

sil<-function(DMAT,labels){
if(nrow(DMAT)!=ncol(DMAT)) DMAT<-eucmetric(t(DMAT))
if(ncol(DMAT) != length(labels)) stop("Matrix dimension and labels length do not match")
sil<-.Call("icav_c_silhouette",DMAT,as.integer(labels),unique(labels),length(labels),DUP=F)
return(sil)
}
