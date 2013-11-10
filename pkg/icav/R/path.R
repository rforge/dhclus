mst2Path.Diss<-function(X){

mst<-function(X){
#X dist matrix or diss matrix --- must be simmetric
	nc<-ncol(X)
	mst<-.Call("icav_mstPath",X,as.integer(nc),DUP=F,PACKAGE="icav")
	return(mst)
}

	MST<-mst(X)
	Diss_1<-.Call("icav_pathDiss_1",MST,as.integer(c(length(MST[[1]]),0)),DUP=F,PACKAGE="icav")

	return(Diss_1)
}
