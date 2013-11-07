Get.clusters<-function(data, index, diss=FALSE, debug=FALSE,method=2,...)
{
  
  if(method==1) {
    s<-select_k(data[index,],kmax=9,Ca=7,centers=2,debug=debug)
    val<- s$k #sel2_k(data[index,])
    
    if(debug)
    {
      print(paste("Selected k: ",val))
    } 
    
    Dist<-as.matrix(pknng(data[index,],k=val,diss=F,fixed.k=1, silent=T, MinGroup=0,penalize=1))
    sigma<-s$sigma 
    psigma<-s$psigma
    Sx1<-exp(-Dist^2/(2*sigma^2))
    if(!s$error){
      sc<-spectral.clust(Sx1, 2) 
      sc$error<-FALSE
      sc$val<-val
      sc$sigma<-sigma
      sc$sigmas<-psigma
    }
    else 
    {
      sc<-list("")
      print("spectral error")
      sc$specc.error<-TRUE
      sc$error<-TRUE
    }
  }
  
  else 
  { print(length(index))
    Dist<-as.matrix(dist(data[index,]))
    
    re<-spectral.sep(Dist,method=method,centers=2,kmax=10,Ca=7,debug=debug)
    sc<-re[[1]]
    Sx1<-re[[2]]
    sigma<-re[[3]]
    psigma<-re[[4]]
    sc$val<-2
    sc$sigma<-sigma
    sc$sigmas<-psigma
    sc$error<-FALSE
  } 
  
  if(!sc$specc.error){
    outliers<-index[sc$is.outlier]
    
    if(length(outliers)>0) { index <-index[!sc$is.outlier]
                             tags[outliers]<-0
                             if (debug){
                               print("Outliers found:")
                               print(outliers)
                             }
                             sc$clusters[sc$is.outlier]<-0                
                             
                             #relabel is needed if outliers are found
                             sc$clusters<-relabel(sc$clusters[!sc$is.outlier])
                             indi<- 1:dim(sc$yi)[1]
                             ind<-indi[!sc$is.outlier] 
                             
                             
    }
  }else
    sc$error<-TRUE
  return(sc)
}



#######################
first_nonzero<-function(v,k){
  flag<-0
  while(v[k]==0){
    k<-k+1
    if(k>=length(v)){ flag<-1;break;}
  }
  if(flag==1) return (1)
  else return (v[k])
}


#### Make Local Scaling similarity matrix
rbf.dot.multiscale2 <-function(Dist,sigma)
{ 
  W<-matrix(0, ncol=dim(Dist)[2], nrow=dim(Dist)[1])
  for(i in 1:dim(Dist)[1])
  { W[i,]<-(-Dist^2)[i,]/(2*sigma*sigma[i])   
  } 
  return(exp(W))
}

#### Make Local Scaling similarity matrix
# rbf.dot.multiscale<-function(data,scale.vector){
#   if(dim(data)[2]!=length(scale.vector)) stop("Length Missmatch")
#   euc<-.Call("rbf_dot_multiscaled",data, c(dim(data)[1],dim(data)[2]),as.double(scale.vector),DUP=F)
#   return(euc)
# }

##Select sigma for gaussian kernel
select.sigma<-function(Dist,method,kmax=10,centers,Ca,debug)
{
  kmax=min(kmax,floor(dim(Dist)[1]/2)+1)
  Ca<-8
  error<-0
  if(method!=3){
    pots<-1/c(sapply(1:Ca,function(x)(2^x)))  
    if (method==1) pots<-c(9:6/10,pots)[1:8]
    print(pots)  
    tmu<-quantile(Dist,probs=pots)
    C<-length(tmu)
  }
  else { 
    sort.X<-apply(Dist,1,sort) 
    if(dim(Dist)[1] > (kmax+2)){ 
      C<-kmax
      tmu<-rep(1,kmax)
      ttmu<-sort.X[3:(kmax+2),]    
      
    }
    else 
    {
      if (dim(Dist)[1] > 3)
      {C<-1
       tmu<-0
       ttmu<-sort.X[2,]}
      else {error <-1 
            print("error")}
    }}
    
    if(!error){
      diss <-(rep(Inf,C))
      w2 <-(rep(Inf,C))
      w1 <-(rep(Inf,C))
      bal<-(rep(0,C))
      i<-1
      
      while (i <C)
      {
        if(method==3) {
          ka <- rbf.dot.multiscale2(Dist,ttmu[i,])
          
        }
        else ka <-exp((-(Dist^2))/(2*(tmu[i]^2)))
        diag(ka) <- 0  
        d <- 1/(rowSums(ka))
        
        if(!any(d==Inf) && !any(is.na(d)))# && (max(d)[1]-min(d)[1] < 10^4))
        {    
          res <- spectral.clust(ka,centers)
          
          if(sum(res$is.outlier)>0){
            res$clusters<-relabel(res$clusters)
            indi<- 1:dim(res$yi)[1]
            ind<-indi[!res$is.outlier]   
            print("outliers:")
            print(indi[res$is.outlier])
            ind<-ind[!is.na(ind)]
            res$yi<-res$yi[ind,]
            res$clusters<-relabel(res$clusters[ind])
            
          }else  ind<- 1:dim(res$yi)[1]
          
          if(length(unique(res$clusters))==centers){
            
            ws2<-sum(withinsum2(res$yi[,1:(centers)],res$clusters))
            ws1<-withinsum2(res$yi[,1:(centers)],rep(1,length(res$clusters)))
            w2[i]<-ws2
            w1[i]<-ws1
            diss[i] <-ws2 /abs(ws1-ws2)
            
            li<-rep(0,centers)
            for(h in 1:centers) li[h]<-sum(res$clusters==h)    
            bal[i]=min(li)/max(li)
            if(min(li)<=2)
            {
              print("singletons detected")
              diss[i]<-Inf
            }
            i<-i+1
          } 
          else 
          {print("Del")
           i<-i+1     
           
          } }
        else 
        {   print(paste ("Deleted2 ",i))
            i<-i+1}
      }
      
      
      ms <- which.min(diss )
      
      pms<-ms
      print(pms)
      if(debug){
        print(diss)
        print(bal)
        if (method!=3)
          print(paste("selected ",ms,"sigma: ",tmu[ms]))
        
      }
#       if(method!=3)
#       {
#         pots2<-c((   c(-1,-1,-1,1,1,1)/(2^c(13,11,9,9,11,13))) + pots[ms],pots[ms]) 
#         print("segundo")
#         print(pots2)
#         tmu<-c(quantile(Dist,probs=pots2))
#         diss <-(rep(Inf,length(tmu)))
#         w2 <-(rep(Inf,length(tmu)))
#         w1 <-(rep(Inf,length(tmu)))
#         bal<-(rep(0,length(tmu)))
#         i<-1
#         #segundo barrido
#         
#         while (i <length(tmu)+1){
#           ka <-exp((-(Dist^2))/(2*(tmu[i]^2)))
#           diag(ka) <- 0  
#           d <- 1/(rowSums(ka))
#           
#           if(!any(d==Inf) && !any(is.na(d))) #&& (max(d)[1]-min(d)[1] < 10^4))
#           {    
#             res <- spectral.clust(ka,centers)
#             
#             if(sum(res$is.outlier)>0){
#               res$clusters<-relabel(res$clusters)
#               indi<- 1:dim(res$yi)[1]
#               ind<-indi[!res$is.outlier]   
#               print("outliers:")
#               print(indi[res$is.outlier])
#               ind<-ind[!is.na(ind)]
#               res$yi<-res$yi[ind,]
#               res$clusters<-relabel(res$clusters[ind])
#               
#             }else  ind<- 1:dim(res$yi)[1]
#             
#             if(length(unique(res$clusters))==centers){
#               
#               ws2<-sum(withinsum2(res$yi[,1:centers],res$clusters))
#               ws1<-withinsum2(res$yi[,1:centers],rep(1,length(res$clusters)))
#               w2[i]<-ws2
#               w1[i]<-ws1
#               diss[i] <-ws2 #/abs(ws1-ws2)
#               li<-rep(0,centers)
#               for(h in 1:centers)
#                 li[h]<-sum(res$clusters==h)
#               bal[i]=min(li)/max(li)
#               if(min(li)<=2)
#               {
#                 print("singleton detected")
#                 diss[i]<-Inf
#               }
#               i<-i+1
#             } 
#             else 
#             {print("Del")
#              i<-i+1     
#              
#             } }
#           else 
#           {   print(paste ("Deleted2 ",i))
#               i<-i+1}
#           
#         }
#         
#       }
#       if(method!=3){
#         ms <- which.min(diss )
#       if(unique(bal)>1) print(":::Cambia") else print(":::NOcambia")
#       }
#       
      if (method==3) {
        if(dim(Dist)[1] > (kmax+2)) mu <-ttmu[(pms),]
        else {mu<-ttmu
              pms<-0}
      }
      else{ ms<-which.min(diss  )
            mu<-tmu[ms]
            if(debug){  
              print(diss)
              print(bal)
              print(paste("selected ",ms," sigma:",tmu[ms]))
              print("......................")
            }
      }
      if(is.infinite(diss[ms])) error<-1
      else error<-0
      
      out<-list(diss,mu,tmu, diss[ms],ms,bal[ms],error,pms) 
      names(out)<-list("diss","mu","tmu", "min","ms","bal","error","pms")
    } else out<-NULL
    return(out)
    
  }
  
##Spectral clustering for method 2 and 3
## Dist is a distance matrix
spectral.sep<-function(Dist,method,centers,kmax,Ca,debug)
{
  s<-select.sigma(Dist, method=method,centers=centers,kmax=kmax,Ca=Ca,debug=debug)
  if(!is.null(s)){
  tmu<-s$tmu
  diss<-s$diss
  
  if(length(diss[diss==Inf]) < length(tmu) && !s$error ) 
  {   
    mu<-s$mu

    if(method==3){
      Sx <- rbf.dot.multiscale2(Dist,mu)   
    }
    else{
      Sx<-exp(-Dist^2/(2*mu^2))
      
    }
    sc<-spectral.clust(Sx, 2)  
    
    psigma<-s$pms
  } 
  else
  {
    sc<-list("")
    print("spectral error")
    sc$yi<-list("")
    sc$specc.error<-TRUE
    Sx<-list("")
    mu<-1
    psigma<-1
    
  }}
  else 
  {
    sc<-list("")
    print("spectral error")
    sc$yi<-list("")
    sc$specc.error<-TRUE
    Sx<-list("")
    mu<-1
    psigma<-1
    
  }
  return(list(sc,Sx,mu,psigma))

    
}


##Select k and sigma for method 1
select_k<-function(data, diss=FALSE,kmax=15,Ca,centers,debug)
{
  kmin<-2
  l<-dim(data)[1]
  kmax<- min(kmax,floor((l)/2) )
  if(kmax > kmin){
    cant<-(kmax-kmin)+1
    sigmas<-vector(length=cant)
    psigmas<-vector(length=cant)
    val<-vector(length=cant)
    bal<-vector(length=cant)
    error<-vector(length=cant)
    i<-1
    while(i<cant+1)
    {
      Dist<-as.matrix(pknng(data,k=((kmin+i)-1),diss=diss,fixed.k=1, silent=T, MinGroup=0))
      print(kmin+i-1)
      s<-select.sigma(Dist,method=1,centers=centers,Ca=Ca,debug=debug)
      sigmas[i]<-s$mu
      psigmas[i]<-s$pms
      val[i]<-s$min
      bal[i]<-s$bal
      error[i]<-s$error
  #    if(bal[i]==1) i<-cant+1
      i<-i+1
    }
    
    rm(Dist) 
    print(paste("dim ",l))
    ks<-which.min(val)
    sigma<-sigmas[ks]
    psigma<-psigmas[ks]
    print(paste("validity:", val))   
    err<-error[ks]
    ks<-ks+1
  }
  else {
    ks<-kmin
    Dist<-as.matrix(pknng(data,k=3,diss=diss,fixed.k=1, silent=T, MinGroup=0))
    s<-select.sigma(Dist,method=1,centers=centers,Ca=Ca,debug=debug)
    sigma<-s$mu
    psigma<-s$pms
    err<-s$error
  }
  
  if(debug){
    print(paste("ks:", ks))
    print(sigma)
  }
  out<-list(ks, sigma,err,psigma)
  names(out)<-list("k","sigma","error","psigma")
  return(out)
}


##Spectral clustering with Lrw Laplacian
## W is a similarity matrix

spectral.clust<-function(W, centers){
  if (length(centers)==1) nc<-centers
  else nc<-dim(centers)[1]
  BIG.NUM<-10e100
  N<-dim(W)[1]
  ERROR.SPECC<-FALSE 
  diag(W)<-0
  d<-1/rowSums(W)
  
  bool.inf<-rep(FALSE,N)
  if(any(d==Inf)  ){
    bool.inf<-is.infinite(d)
    
    if(any(bool.inf)){
      inf.clusters<-rep(0,N)
      inf.clusters[bool.inf]<-1:sum(bool.inf)
      W<-W[!bool.inf,!bool.inf]
      d<-d[!bool.inf]
    }   
  }
  
  
  if(!any(d==Inf) && !any(is.na(d)) && length(d)>0)
  {   
    L.rw<- diag(d)%*%W 
    ss<-eigen(L.rw,symmetric=F)$vectors[,1:nc]   
    yi <- Re(ss)
    rm(W,ss)
  }
  else{ERROR.SPECC<-TRUE;yi<-0}
  
  if (any(is.na(yi)) || any(is.infinite(yi)) || ERROR.SPECC){
    res<-list()
    print("Spectral Error")
    res$specc.error<-TRUE
    res$singletons<-FALSE
    #    res$sum.specc.withinss<- Inf
    res$cluster<-rep(1,N)
    res$clusters<-res$cluster
    res$yi<-yi
    res$is.outlier<-bool.inf
    res$affinity.matrix<-L.rw
  }
  else{
    res<- pam(dist(yi[,1:nc]),nc)
    #    res$sum.specc.withinss<- sum(res$clusinfo[,3])
    clusters<-res$clustering    
    res$singletons<-FALSE
    
    if(any(bool.inf)){
      inf.clusters[!bool.inf]<-res$clustering+sum(bool.inf)
      clusters<-inf.clusters
      res$singletons<-TRUE
    }    
    res$is.outlier<-bool.inf
    res$clusters<-clusters
    res$specc.error<-FALSE
    res$yi<-yi
    res$affinity.matrix<-L.rw
  }
  return(res)
}

