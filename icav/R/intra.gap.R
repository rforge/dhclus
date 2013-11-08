# icav.estimator.apply<-function(labels=NULL, x,filter.outliers=0,dump=F){
# 	if (is.null(labels)) stop("Check labels")
# 		
# 		sort.D<-function(data){return(sort(data,index.return = TRUE)$ix)}
# 		quant.D<-function(data){return(quantile(data,probs = c(0.25,0.5,0.75) ))}
# 		D<-as.matrix(dist(x))
# 		diag(D)<-Inf
# 		index<-1:nrow(x)
# 		
# 		
# 		DMIN<-list()
# 		DMIN$dist<-0
# 		DMIN$i<- -1
# 		DMIN$j<- -1
# 		DMIN$outliers<-list()
# 		index.list<-list()
# 		u.labels<-unique(labels)
# 		intra.class<-vector()
# 		
# 		
# 		if (filter.outliers==0) intra.knn<-1
# 			else intra.knn<-filter.outliers
# 				for(i in 1:length(u.labels)){
# 					index.1<-index[labels==u.labels[i]]
# 					
# 					if(filter.outliers>0 && length(index.1) > 9){## optimizar para hacerlo una vez
# 						D.intra.1<-icav.Matrix.chunk(D,index.1,index.1)
# 						D.sort.i1<-icav.Matrix.sort(t(D.intra.1))
# 						D.stat.i1<-quant.D(D.sort.i1[filter.outliers,])
# 						out.12<-D.stat.i1[2]+(D.stat.i1[3]-D.stat.i1[1])*1.5 # 0.75 quantile
# 						bool.1<-D.sort.i1[filter.outliers,]< out.12
# 						DMIN$outliers[[i]]<-index.1[!bool.1]
# 						index.list[[i]]<-index.1[bool.1]
# 						
# 						
# 					}
# 					else index.list[[i]]<-index.1
# 						
# 						if (length(index.list[[i]])>1) {
# 							D.intra.1<-mst2Path.Diss(icav.Matrix.chunk(D,index.1,index.1)); 
# 							DMIN$dist<-DMIN$dist + .Call("icav_wss_matrix",D.intra.1,dim(D.intra.1),DUP=F)
# 						}
# 				}
# 				
# 				return(DMIN$dist)
# }
# icav.estimator<-function(x, labels=NULL, dump=F){
# if (is.null(labels)) stop("Check labels")
# 
# sort.D<-function(data){return(sort(data,index.return = TRUE)$ix)}
# quant.D<-function(data){return(quantile(data,probs = c(0.25,0.5,0.75) ))}
# D<-as.matrix(dist(x))
# diag(D)<-Inf
# index<-1:nrow(x)
# 
# 
# DMIN<-list()
# DMIN$dist<-0
# DMIN$i<- -1
# DMIN$j<- -1
# DMIN$outliers<-list()
# index.list<-list()
# u.labels<-unique(labels)
# intra.class<-vector()
# 
# 
# if (filter.outliers==0) intra.knn<-1
# else intra.knn<-filter.outliers
# for(i in 1:length(u.labels)){
#     index.1<-index[labels==u.labels[i]]
# 
#     if(filter.outliers>0 && length(index.1) > 9){## optimizar para hacerlo una vez
#       D.intra.1<-icav.Matrix.chunk(D,index.1,index.1)
#       D.sort.i1<-icav.Matrix.sort(t(D.intra.1))
#       D.stat.i1<-quant.D(D.sort.i1[filter.outliers,])
#       out.12<-D.stat.i1[2]+(D.stat.i1[3]-D.stat.i1[1])*1.5 # 0.75 quantile
#       bool.1<-D.sort.i1[filter.outliers,]< out.12
#       DMIN$outliers[[i]]<-index.1[!bool.1]
#       index.list[[i]]<-index.1[bool.1]
# 
# 
#     }
#     else index.list[[i]]<-index.1
# 
#     if (length(index.list[[i]])>1) {
# 	D.intra.1<-mst2Path.Diss(icav.Matrix.chunk(D,index.1,index.1)); 
# 	DMIN$dist<-DMIN$dist + .Call("icav_wss_matrix",D.intra.1,dim(D.intra.1),DUP=F)
#     }
# }
# 
#   return(DMIN)
# }


icav.estimator.apply<-function(labels=NULL, x,dump=F){
	if (is.null(labels)) stop("Check labels")
	DMIN<-0
# 	print(labels)
	D<-as.matrix(dist(x))
	D.intra.1<-mst2Path.Diss(D); 
	DMIN<-sum(.Call("icav_gap_withinSum",D.intra.1,as.integer(labels),unique(labels),length(labels),DUP=F))

	return(DMIN)
}

icav.estimator<-function(x, labels=NULL, dump=F){
if (is.null(labels)) stop("Check labels")

	D<-as.matrix(dist(x))
	D.intra.1<-mst2Path.Diss(D); 
	DMIN<-sum(.Call("icav_gap_withinSum",D.intra.1,as.integer(labels),unique(labels),length(labels),DUP=F))
    


  return(DMIN)
}


icav.gap<-function(data, clustering=NULL, k.max = 20, M = 100){

x<-data
k.max <- trunc(k.max)
if(k.max < 2) stop("'k.max' has to be >= 2")

    clust.obs<- clustering(x, k.max = k.max)

    D.obs<-vector()
    D.obs[1]<-icav.estimator.apply(labels=rep(1,nrow(x)),x )
    for (i in 1:(k.max-1)) {D.obs[i+1]<-icav.estimator.apply(labels=clust.obs[,i],x )}

  x<-scale(x,center=TRUE,scale=FALSE)
  veigen <- svd(x)$v
  x1 <- crossprod(t(x), veigen)
  
  N<-ceiling(nrow(x1))
  z1 <- matrix(data = NA, nrow = N, ncol = ncol(x1))

  min.x <- apply(x1,2,min)
  max.x <- apply(x1,2,max)
  
  clust.null<-list()
D.null<-matrix(0,M,k.max)
  tots <- vector(length = M)
  for (k in 1:M) {
    for (j in 1:ncol(x1)) {
      z1[, j] <- runif(N, min = min.x[j], max = max.x[j])
    }
	  clust.null[[k]]<-list()
	  cl<-clustering(z1, k.max = k.max)

	D.null[k,1]<-icav.estimator.apply(labels=rep(1,nrow(z1)), z1 )

	  for (i in 1:(k.max-1)) {
	    D.null[k,i+1]<-icav.estimator(labels=cl[,i], z1)
	  }
  }
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



res<-list()
res$D.obs<-D.obs
res$clust.obs<-clust.obs
res$D.null<-D.null
gapStat<-cbind(1:length(D.obs),apply(log10(D.null),2,mean)-log10(D.obs),sqrt(1 + 1/M) * apply(log10(D.null),2,sd),apply(D.null,2,mean), D.obs)
colnames(gapStat) <- c("No. of Clusters","Gap statistic", "SE of simulation", "Null Withiness", "Obs. Withiness")
res$gapStat<-gapStat

return(res)
}

icav.gap.mc<-function(data, clustering=NULL, k.max = 20, M = 100, filter.outliers=0,mc.cores=8){

x<-data
k.max <- trunc(k.max)
if(k.max < 2) stop("'k.max' has to be >= 2")

    clust.obs<- clustering(x, k.max = k.max)

    D.obs<-vector()
    D.obs[1]<-icav.estimator.apply(labels=rep(1,nrow(x)),x )
    for (i in 1:(k.max-1)) {D.obs[i+1]<-icav.estimator.apply(labels=clust.obs[,i],x )}
  
  x<-scale(x,center=TRUE,scale=FALSE)
  veigen <- svd(x)$v
  x1 <- crossprod(t(x), veigen)
  
  N<-ceiling(nrow(x1))

  min.x <- apply(x1,2,min)
  max.x <- apply(x1,2,max)
  clust.null<-list()
  tots <- vector(length = M)
  
  do.rand.matrix<-function(dim.i,N.data,min.val,max.val){return(runif(N.data,min.val[dim.i],max.val[dim.i]))}
  
  foo<-function(val,N,min.x,max.x,kmax) {
	z1<- matrix( unlist(lapply(1:length(min.x), do.rand.matrix, N, min.x, max.x) ) , nc=length(min.x), nr=N )
	cl<-clustering(z1, k.max = k.max)
	cl.val<-apply(cbind(1,cl), 2, icav.estimator.apply, z1)
	return(cl.val)
  }
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
clust.null<-mclapply(1:M,foo,N,min.x,max.x,kmax,mc.cores=mc.cores)

D.null<- t(matrix(unlist(clust.null),nc=M,nr=k.max))
    
res<-list()

res$D.obs<-D.obs
res$clust.obs<-clust.obs
res$D.null<-D.null
res$D.1knn.null<-NULL
res$D.1knn<-NULL
gapStat<-cbind(1:length(D.obs),apply(log10(D.null),2,mean)-log10(D.obs),sqrt(1 + 1/M) * apply(log10(D.null),2,sd),apply(D.null,2,mean), D.obs)
colnames(gapStat) <- c("No. of Clusters","Gap statistic", "SE of simulation", "Null Withiness", "Obs. Withiness")
res$gapStat<-gapStat
return(res)
}

