pknng<-function(X,k,fixed.k=0,diss=F,cte=3,mu="mean",method="euclidean",conn="one",penalize=1,path="dijkstra",MinGroup=0,silent=T,SMP=T){
# penalize=1 exp, penalize=2 lineal ,penalize=3 power

data<-t(X)
ret<-list()
if (diss) X<-as.matrix(X)
else{
if (method=="euclidean") {
if (SMP) X<-eucmetric(data)
else X<-as.matrix(dist(t(data)))
}else if (method=="corr") {

X<-2*cormetric(data)
X[X<0]<-0
X<-sqrt(2*X)
}
else if (method=="mi.1") X<-mi_metric(data)
else if (method=="mi.2") X<-mi_estimator(data)
else stop("method options: euclidean,corr,mi.1,mi.2")
}
rm(data)
# print (paste("Data row: ",dim(X)[1],"by col: ",dim(X)[2]))
m <- nrow(X)
if(k>m){stop("k > number of observations")}
if(!silent) print (paste("Data row: ",dim(X)[1],"by col: ",dim(X)[2]))

thr<- -1

if (conn=="one") conntype<-1
else if (conn=="ttr") conntype<-2
else if (conn=="all") conntype<-0
else stop("Error\n")
 if (conn=="all" && penalize==0) stop("connection type all cannot have penalization 0 (non-penalized all matrix is a distance matrix)\n")

# print("get_k_neighbors")

ret<-.Call("get_k_neighbors",X,c(dim(X)[1],dim(X)[2],k),DUP=F)

# print("fin get_k_neighbors")
# print(dim(X))
len.ret<-length(ret[[1]])
v.ret<-unlist(ret[[1]])
## no borrar se usa abajo mu = 1q y mu = 3q
qq<-quantile(v.ret,probs=c(0.25,0.75)) 
rm(v.ret)
sg<-qq[2]+(qq[2]-qq[1])*1.5

ret<-.Call("make_symmetric",as.list(ret),c(dim(X)[1],k),sg,DUP=F)

#  cat("\n\n")
#  print(qq)
#  print(sg)
#  print(ret)

# print(dim(ret[[1]]))
#print(dim(X))
#  print("make_symmetric")
if (fixed.k==0){
#  print("make_symmetric")
v<-unlist(ret[[1]])
if (mu=="min") mu<-min(v)
if (mu=="median") mu<-median(v)
else if (mu=="mean") mu<-mean(v)
else if (mu=="max") mu<-max(v)
else if (mu=="1q") mu<-qq[1]
else if (mu=="3q") mu<-qq[2]
rm(v)
}
else{
if(!silent)  print(dim(X))
sort.X<-apply(X,1,sort)
mu<-apply(sort.X,1,quantile,probs=c(0.75))
mu<-mu[fixed.k+1]
# print("mu: ")
# print(mu)
}
if(!silent) {print("mu: ");print(mu)}

if (!is.list(ret)) stop("Data Problem\n")
metric<-0

if (MinGroup > 0){
addons<-.Call("outlayer_merger",as.list(ret),X,as.integer(c(dim(X)[1], dim(X)[2], conntype, penalize, metric, cte, MinGroup)),DUP=F)
len<-length(addons[[length(addons)]])
if (len>1){
conns<-as.matrix(addons[[length(addons)]])
for(i in 1:dim(conns)[1]){
ret[[2]][[conns[i,1]]]<-as.integer(c(ret[[2]][[conns[i,1]]],conns[i,2]))
ret[[2]][[conns[i,2]]]<-as.integer(c(ret[[2]][[conns[i,2]]],conns[i,1]))
ret[[1]][[conns[i,1]]]<-as.double(c(ret[[1]][[conns[i,1]]],conns[i,3]))
ret[[1]][[conns[i,2]]]<-as.double(c(ret[[1]][[conns[i,2]]],conns[i,3]))
}
rm(addons,len)
}
}
# print(paste("mu: ",mu))
# cat("\n\n")
# print(sg)
# print("connect")
MinGroup<-0
if (conn!="all"){
if (conn=="one") addons<-.Call("connect_groups",as.list(ret),X,as.integer(c(dim(X)[1], dim(X)[2], conntype, penalize, metric, cte, MinGroup )),as.double(mu),DUP=F)
else if (conn=="ttr"){
 addons<-.Call("connect_groups_oneToTheRest",as.list(ret),X,as.integer(c(dim(X)[1],dim(X)[2],conntype,penalize,metric,cte)),as.double(mu),DUP=F)
}
# print("fin connect")
# print(addons)
len<-length(addons[[length(addons)]])
# print(len)

if (len>1){
conns<-as.matrix(addons[[length(addons)]])
# print(addons)
for(i in 1:dim(conns)[1]){
ret[[2]][[conns[i,1]]]<-as.integer(c(ret[[2]][[conns[i,1]]],conns[i,2]))
ret[[2]][[conns[i,2]]]<-as.integer(c(ret[[2]][[conns[i,2]]],conns[i,1]))
ret[[1]][[conns[i,1]]]<-as.double(c(ret[[1]][[conns[i,1]]],conns[i,3]))
ret[[1]][[conns[i,2]]]<-as.double(c(ret[[1]][[conns[i,2]]],conns[i,3]))
}
}
}
else addons<-.Call("connect_groups_all",as.list(ret),X,as.integer(c(dim(X)[1],dim(X)[2],conntype,penalize,metric,cte)),as.double(mu),DUP=F)
rm(X)
# print(ret[[2]])
# cat("\n\n")
# print("all_dijkstra")
if(path=="dijkstra")d_matrix<-.Call("all_dijkstra",as.list(ret),c(m),DUP=F)
else if(path=="floyd"){
d_matrix<-.Call("pknng_floyd",as.list(ret),as.integer(c(m,0)),DUP=F)
}
else stop("error\n path == dijkstra or floyd")
# print("fin all_dijkstra")
# print(d_matrix)
return(d_matrix)
}
