
###Test functions


test<-function(data,result,...)
{
  return(test_gap(data,labels=result$clusters,sigma=result$sigma,sigmas=result$sigmas,k=result$val,...))
}


#test gap statistics 
# metric is used for compute the WSS : 1 pknng, 2 icav, 3 euclidean

test_gap<- function(data, labels, NumRef=150,sigmas,k=3,sigma=1,  method=2, debug=T,metric=2)
{ 

 debug=TRUE
  B<-NumRef
  W<-vector()
  WSS<-vector() 
  dat <- scale(data, center = TRUE, scale = F)
 veigen <- svd(dat)$v
 data <-crossprod(t(dat), (veigen))
 x1<-data
 N<-nrow(x1) 
 dat<-x1
  val<-(dim(data)[1])/2
#  if(method==1) s<-k
#  else 
    s<-min(k,val)

    dis <-if(metric==1) pknng(dat,k=s,diss=F,fixed.k=1, silent=T, MinGroup=0,penalize=1)
           else if(metric==2) mst2Path.Diss(as.matrix(dist(dat)))
              else as.matrix(dist(dat))
  ##test
##dis<- as.matrix(dist(dat))  
  D<-dim(dis)[1]

  WSS[1]<- (withinsum(dis, rep(1,D), diss=T))
  W[1]<- log(WSS[1])  
  WSS[2]<-sum(withinsum(dis, labels, diss=T))
  W[2]<- log(WSS[2])
  
 

  min.x <- apply(x1,2,min)
  max.x<-apply(x1,2,max)
  
  foo<-function(v,N,min.x,max.x,veigen,method,sigmas,k,sigma,metric){
    Ws<-vector()
    z11 <- matrix(data = 0, nrow = N, ncol = length(min.x))  
    
    for (j in 1:length(min.x)) 
    {
      z11[,j] <- runif(N, min = min.x[j], max = max.x[j])
    }
    
    X1<-z11##%*%t(veigen)  # back transformed
    dx<- dist(X1)
    dx1<-as.matrix(dx)
    val<-(dim(dx1)[1])/2
    #if(method==1) s<-k
    #else
      s<-min(k,val)
    
    dx<-if(metric==1) pknng(X1,k=s,diss=F,fixed.k=1, silent=T, MinGroup=0,penalize=1)
      else if(metric==2) mst2Path.Diss(as.matrix(dist(X1)))
       else as.matrix(dist(X1))
    
    
    if (method!=3){
     
      if(method==1 ) dx1<-pknng(X1,k=k,diss=F,fixed.k=1, silent=T, MinGroup=0,penalize=1)
      
      pots<-1/c(sapply(1:8,function(x)(2^x)))  
      pots<-(c(9:6/10,pots))#[1:8]
     # if (method==1) pots<-pots[1:7]
      
      ms<-sigmas[[1]]
      sel<-pots[ms]
     # print(sel)
      if(ms==1) a<-(1-sel)/5
      else a<-(pots[(ms-1)]-sel)/(ms+4)
      if(ms==length(pots)) b<-sel/2
      else b<-(sel - pots[(ms+1)])/(ms+2)
      
      pots2<-c(sel+2*a, sel+a,sel,sel-b,sel-2*b) 
      
  
      tmu<-quantile(dx1,probs=pots2)[sigmas[[2]]]
      Sx1<-exp(-dx1^2/(2*tmu^2))
      #   Sx1<-exp(-dx^2/(2*sigma^2))
    }
    else 
    { 
      sort.X<-apply(dx1,1,sort) 
      vec<-min(7,(dim(dx1)[1])/2)
      ttmu<-sort.X[(sigmas+2),]  
      Sx1 <- rbf.dot.multiscale2(dx1,ttmu)
    #  Sx1 <- rbf.dot.multiscale2(dx,sigmas)
    }
    if(method<4){
    labs<-spectral.clust(Sx1, 2)
    if(!labs$specc.error){
    iind<-1:length(labs$yi[,2])
    outliers<-iind[labs$is.outlier]}
    }
    else {
      labs<-list()
      labs$clusters<-kmeans(X1,2,iter.max=200)$cluster
      outliers<-list()
      labs$specc.error<-FALSE
    }
    if(!labs$specc.error&&length(outliers)>0) 
    { Ws<-NA
    } 
    else 
    { 
    Ws[1]<- log(withinsum(dx, rep(1,ncol(dx)), diss=T))
    Ws[2]<- log(sum(withinsum(dx, labs$clusters,diss=T)))
   # print(sum(labs$clusters==1))
    }
    rm(X1,z11,dx)
    # gc()
    return(as.double(Ws))
  }
  
  wss.null<-matrix(0, nrow = NumRef, ncol = 2)
  
  if(co<-detectCores() ){  
   if(debug) print(paste("cores ",co))
    if (Sys.info()[1] == "Windows"){
      cl <- makeCluster(cores)
      l.gap <-clusterApply(cl=cl,1:NumRef,foo,N=N,min.x=min.x,max.x=max.x,veigen=veigen,method=method,
                           sigmas=sigmas,k=k,sigma=sigma, metric=metric)
    }
    else  {l.gap<-multicore::mclapply(1:NumRef,foo,N=N,min.x=min.x,max.x,method=method, sigmas=sigmas,k=k,
                                      sigma=sigma,metric=metric,mc.cores = detectCores() )
           }  
    for (b in 1:NumRef){
      wss.null[b,]<-as.double(unlist(l.gap[[b]]))
    } 
  }
  else{
    for (i in 1:NumRef){
      wss.null[i,]<-foo(i,N,min.x,max.x,veigen=veigen,method=method, sigmas=sigmas,k=k,sigma=sigma, metric=metric)   
    }
}
  
  Elogw<-apply(wss.null,2,mean,na.rm=TRUE)
  GAP<-Elogw-W 
  GAP.sd <- sqrt(1 + 1/NumRef) * apply(wss.null,2,sd,na.rm=TRUE)
  
  cond<-GAP[1] <GAP[2]-GAP.sd[2]
  gap1<-abs(GAP[1]  -( GAP[2] -GAP.sd[2]))
  
  out <- list(1:2,W, Elogw, GAP, GAP.sd ,cond,gap1)
  names(out)<-list("Cl","OlogW","ElogW", "GAP", "GAP.sd", "cond", "gap")
  if(debug) print((out))
  return (out)
}




