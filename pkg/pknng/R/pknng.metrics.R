cormetric<-function(data){

correl<-.Call("correlation_metric",data, c(dim(data)[1],dim(data)[2]),DUP=F,PACKAGE="pknng")
#par[0] rows (dimension of Data), par[1] cols (amount of Data)
correl<-matrix(correl,dim(data)[2],dim(data)[2])
return(correl)
}

eucmetric<-function(data){

euc<-.Call("euc_metric",data, c(dim(data)[1],dim(data)[2]),DUP=F,PACKAGE="pknng")
#par[0] rows (dimension of Data), par[1] cols (amount of Data)
return(euc)
}

histfun<-function(data){

euc<-.Call("histo_fun",data,length(data),DUP=F,PACKAGE="pknng")
#par[0] rows (dimension of Data), par[1] cols (amount of Data)
return(euc)
}

mi_metric<-function(data){
euc<-.Call("mi_metric",data,c(dim(data)[2],dim(data)[1],0,1),DUP=F,PACKAGE="pknng")
return(euc)
}

mutualInf<-function(data){
euc<-.Call("mutualInfo",data,c(dim(data)[2],dim(data)[1],0),DUP=F,PACKAGE="pknng")
return(euc)
}

mi_estimator<-function(data){
## data rows dimension, col amount of data
euc<-.Call("mi_estimator",data,c(dim(data)[1],dim(data)[2]),DUP=F,PACKAGE="pknng")
return(euc)
}

rbf.dot.multiscale<-function(data,scale.vector){
if(dim(data)[2]!=length(scale.vector)) stop("Length Missmatch")
euc<-.Call("rbf_dot_multiscaled",data, c(dim(data)[1],dim(data)[2]),as.double(scale.vector),DUP=F,PACKAGE="pknng")
return(euc)
}

multiscale.rbf<-function(x, kpar = 3){
X<-eucmetric(t(x))
ret<-.Call("get_k_neighbors",X,c(dim(X)[1],dim(X)[2],kpar),DUP=F,PACKAGE="pknng")
sigma<-unlist(lapply(ret[[1]],max))

W<-rbf.dot.multiscale(t(x),sigma)
return(W)
}
