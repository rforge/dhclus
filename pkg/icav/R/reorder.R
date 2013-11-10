## c code reorder.c
icav.Matrix.chunk<-function(M,row,col){
size.M<-dim(M)
m<-.Call("icav_Matrix_chunk",M,as.integer(col),as.integer(row),c(size.M, floor(size.M[2]/4) ),DUP=F,PACKAGE=icav)
return(m)
}

icav.Matrix.sort<-function(M){
size.M<-dim(M)
m<-.Call("icav_Matrix_sort",M,c(size.M, floor(size.M[2]/4) ),DUP=F,PACKAGE=icav)
return(m)
}
