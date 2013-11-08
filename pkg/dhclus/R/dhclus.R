dhclus<- function(data, diss=FALSE,debug=FALSE,FUNCluster=Get.clusters, FUNTest=test,...)
{
  ind <- 1:dim(data)[1]
  init_clusters <- ind
  init_clusters[ind]<-0
  
  tags<-calclus.spec(data=data, index=ind,  kl=1 , tags=init_clusters, FUNCluster=FUNCluster, FUNTest=FUNTest,...)
  class<-list()
  class$tags<-tags
  vec<-as.matrix(class$tags)
  max<-(max_level(class$tags))
  if(max >0){
    i<-max
    t<-class$tags
    while(i > 1) 
    { 
      out <- labels_level(t,i)
      t<-out[,1]
      i<-i-1
      print(dim(as.matrix(out)))
      print(dim(vec))
      vec<-cbind(out,vec)
      
    }

    r<-dim(vec)[2]
    class$labelsTree<-vec
#     class$as<-sapply(1:r, function(x) (mean(silhouette(class$labelsTree[,x],dmatrix=as.matrix(mst2Path.Diss(as.matrix(dist(data)))))[,3]) ))
#     class$ssb<-sapply(1:r, function(x) (  bss(mst2Path.Diss(as.matrix(dist(data))),class$labelsTree[,x],diss=TRUE )))
#     class$tags_sil<-relabel(class$labelsTree[,which.max(class$as)])
#     class$tags_bss<-relabel(class$labelsTree[,which.min(class$ssb)])
    #class$tags<-relabel(class$tags_bss)
  }
#  else
#   {
#     class$tags_sil<-relabel(class$tags)
#     class$tags_bss<-relabel(class$tags)
#   }
    class$ltags<-class$tags
    class$tags<-relabel(class$tags)
  class$height<-max
  return(class)
}

calclus.spec <- function(data, index, kl , tags,   diss=FALSE, debug=FALSE, FUNCluster, FUNTest,...)
  #1=pknng, 2 = rbf, 3= local scaling,
{
  tags[index]<-kl
  
  sc<-FUNCluster(data,index, ...)
  
  if(!sc$error){
    left <- index[sc$clusters==1]  
    right<- index[sc$clusters==2]
    ind<-sort(c(left,right))
    out<- index[sc$clusters==0]
    tags[left] <- 2*kl 
    tags[right] <- 2*kl+1
    tags[out] <- 0
    error<-0
  }
  else #spectral error
  {
    print("Error")  
    error <-1
  }
  
  
  if(!error && (g<-FUNTest(data[ind,],sc,...))$cond ) 
  { print(length(left))
    
    print("true")
    kll<-kl
    klr<-kl
    if (length(left)>(3))
    {
      bl <- calclus.spec (data,  index=left, 2*kll , tags,  diss=diss,debug=debug,FUNCluster,FUNTest,...) 
      final <-bl
      tags[left]=bl[left]
      rm(bl)
    }else
    { 
      tags[left]<-kl      
    }
    if (length(right)>(3))
    {      
      br <- calclus.spec (data,  index=right,  (2*klr+1) ,tags,diss=diss, debug=debug,FUNCluster,FUNTest,...)
      tags[right]=br[right]
      rm(br)
    }
    else
    { 
      tags[right]<-kl      
    }    
    rm(left,right)
  }
  
  
  else #No more clusters at this level
  {
    if(debug) print("false")
    tags[index]<-kl
    rm(sc)
  }
  
  
  return(tags)
}









