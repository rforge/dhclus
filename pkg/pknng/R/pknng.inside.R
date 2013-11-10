pknng.insider<-function(X,k,fixed.k=0,diss=T,cte=3,mu="mean",method="euclidean",conn="one",penalize=1,MinGroup=0){
# penalize=1 exp, penalize=2 lineal ,penalize=3 power

ret<-list()
if (diss) X<-as.matrix(X)
else{
if (method=="euclidean") X<-eucmetric(t(X))
else if (method=="corr") X<-cormetric(t(X))
else if (method=="mi.1") X<-mi_metric(t(X))
else if (method=="mi.2") X<-mi_estimator(t(X))
else stop("method options: euclidean,corr,mi.1,mi.2")
}
# print (paste("Data row: ",dim(X)[1],"by col: ",dim(X)[2]))
m <- nrow(X)
if(k>m){stop("k > number of observations")}


thr<- -1

if (conn=="one") conntype<-1
else if (conn=="ttr") conntype<-2
else if (conn=="all") conntype<-0
else stop("Error\n")
 if (conn=="all" && penalize==0) stop("connection type all cannot have penalization 0 (non-penalized all matrix is a distance matrix)\n")

# print("get_k_neighbors")

ret<-.Call("get_k_neighbors",X,c(dim(X)[1],dim(X)[2],k),DUP=F,PACKAGE="pknng")

# print("fin get_k_neighbors")
# print(dim(X))

len.ret<-length(ret[[1]])
v.ret<-unlist(ret[[1]])
## no borrar se usa abajo mu = 1q y mu = 3q
qq<-quantile(v.ret,probs=c(0.25,0.75)) 
rm(v.ret)
sg<-qq[2]+(qq[2]-qq[1])*1.5
ret.2<-.Call("make_symmetric",as.list(ret),c(dim(X)[1],k),sg,DUP=F,PACKAGE="pknng")

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
print(dim(X))
sort.X<-apply(X,1,sort)
mu<-apply(sort.X,1,quantile,probs=c(0.75))
mu<-mu[fixed.k+1]
# print("mu: ")
# print(mu)
}


metric<-0
if (MinGroup > 0){
addons<-.Call("outlayer_merger",as.list(ret.2),X,as.integer(c(dim(X)[1], dim(X)[2], conntype, penalize, metric, cte, MinGroup)),DUP=F,PACKAGE="pknng")
len<-length(addons[[length(addons)]])
if (len>1){
conns<-as.matrix(addons[[length(addons)]])
for(i in 1:dim(conns)[1]){
ret.2[[2]][[conns[i,1]]]<-as.integer(c(ret.2[[2]][[conns[i,1]]],conns[i,2]))
ret.2[[2]][[conns[i,2]]]<-as.integer(c(ret.2[[2]][[conns[i,2]]],conns[i,1]))
ret.2[[1]][[conns[i,1]]]<-as.double(c(ret.2[[1]][[conns[i,1]]],conns[i,3]))
ret.2[[1]][[conns[i,2]]]<-as.double(c(ret.2[[1]][[conns[i,2]]],conns[i,3]))
}
# rm(addons,len)
}
}

MinGroup<-0
if (conn!="all"){
if (conn=="one") addons<-.Call("connect_groups",as.list(ret.2),X,as.integer(c(dim(X)[1], dim(X)[2], conntype, penalize, metric, cte, MinGroup )),as.double(mu),DUP=F,PACKAGE="pknng")
else if (conn=="ttr"){
 addons<-.Call("connect_groups_oneToTheRest",as.list(ret.2),X,as.integer(c(dim(X)[1],dim(X)[2],conntype,penalize,metric,cte)),as.double(mu),DUP=F,PACKAGE="pknng")
}
}
else addons<-.Call("connect_groups_all",as.list(ret.2),X,as.integer(c(dim(X)[1],dim(X)[2],conntype,penalize,metric,cte)),as.double(mu),DUP=F,PACKAGE="pknng")


ret<-.Call("connect_kneigbours",as.list(ret.2),as.integer(c(dim(X)[1],dim(X)[2])),DUP=F,PACKAGE="pknng")
c.2<-.Call("tagger",as.list(ret),c(dim(X)[1],0),DUP=F,PACKAGE="pknng")

ret.neig<-list()
ret.neig$tag<-c.2
ret.neig$groups<-ret
ret.neig$adj.list<-ret.2
ret.neig$connections<-addons[[2]]
return(ret.neig)
}
 
pknng.insider.1<-function(X,k,diss=T,cte=3,mu="mean",method="euclidean",conn="one",penalize=1){
# penalize=1 exp, penalize=2 lineal ,penalize=3 power

ret<-list()
if (diss) X<-as.matrix(X)
else{
if (method=="euclidean") X<-eucmetric(t(X))
else if (method=="corr") X<-cormetric(t(X))
else if (method=="mi.1") X<-mi_metric(t(X))
else if (method=="mi.2") X<-mi_estimator(t(X))
else stop("method options: euclidean,corr,mi.1,mi.2")
}
# print (paste("Data row: ",dim(X)[1],"by col: ",dim(X)[2]))
m <- nrow(X)
if(k>m){stop("k > number of observations")}


thr<- -1

if (conn=="one") conntype<-1
else if (conn=="ttr") conntype<-2
else if (conn=="all") conntype<-0
else stop("Error\n")
 if (conn=="all" && penalize==0) stop("connection type all cannot have penalization 0 (non-penalized all matrix is a distance matrix)\n")

# print("get_k_neighbors")

ret<-.Call("get_k_neighbors",X,c(dim(X)[1],dim(X)[2],k),DUP=F,PACKAGE="pknng")

# print("fin get_k_neighbors")
# print(dim(X))
len.ret<-length(ret[[1]])
v.ret<-ret[[1]][[1]]
for(i in 2:len.ret) v.ret<-c(v.ret,ret[[1]][[i]])
## no borrar se usa abajo mu = 1q y mu = 3q
qq<-quantile(v.ret,probs=c(0.25,0.75)) 
rm(v.ret)
sg<-qq[2]+(qq[2]-qq[1])*1.5
#  cat("\n\n")
#  print(qq)
#  print(sg)
#  print(ret)

# print(dim(ret[[1]]))
#print(dim(X))
#  print("make_symmetric")
ret.2<-.Call("make_symmetric",as.list(ret),c(dim(X)[1],k),sg,DUP=F,PACKAGE="pknng")
# print(ret.2[[2]])
# cat("\n\n")
# print(sort(ret.2[[2]][[1]]))
# cat("\n")
# print(min(ret.2[[2]][[1]]))
# print(sort(ret.2[[2]][[min(ret.2[[2]][[1]])]]))
ret<-.Call("connect_kneigbours",as.list(ret.2),as.integer(c(dim(X)[1],dim(X)[2])),DUP=F,PACKAGE="pknng")
c.2<-.Call("tagger",as.list(ret),c(dim(X)[1],0),DUP=F,PACKAGE="pknng")

return(c.2)
}

NVI.list<-function(class.list,total.points,OMP=T){
l<-length(class.list)

if (OMP) VI.list<-.Call("NMI_intersect_l_OMP",as.list(class.list),c(total.points,length(class.list)),DUP=F,PACKAGE="pknng")
else VI.list<-.Call("NMI_intersect_l",as.list(class.list),c(total.points,length(class.list)),DUP=F,PACKAGE="pknng")

return(VI.list)
}
