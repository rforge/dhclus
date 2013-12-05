

dhclus<- function(data, diss=FALSE,debug=FALSE,file="Smatrix.pdf",FUNCluster=Get.clusters, FUNTest=test,...)
{
#  library(lattice)
 
#  trellis.device(dev=pdf, file=file) 
  n<-dim(data)[1]
  ind <- 1:n
  init_clusters <- ind
  init_clusters[ind]<-0
  Sx<-matrix(0,ncol=n, nrow=n)
  res<-calclus.spec(data=data, index=ind,  kl=1 ,Sx, tags=init_clusters, FUNCluster=FUNCluster, FUNTest=FUNTest,...)
  class<-list()
  class$tags<-res$tags
  class$Sx<-res$Sx
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
 #   class$as<-sapply(1:r, function(x) (mean(silhouette(class$labelsTree[,x],dmatrix=as.matrix(mst2Path.Diss(as.matrix(dist(data)))))[,3]) ))
#    class$ssb<-sapply(1:r, function(x) (  bss(mst2Path.Diss(as.matrix(dist(data))),class$labelsTree[,x],diss=TRUE )))
 #   class$tags_sil<-relabel(class$labelsTree[,which.max(class$as)])
#    class$tags_bss<-relabel(class$labelsTree[,which.min(class$ssb)])
    #class$tags<-relabel(class$tags_bss)
  }
    class$ltags<-class$tags
    class$tags<-relabel(class$tags)
  class$height<-max
 # pdf(file=file,paper="special",width=10,height=10)
  #trellis.device(dev=pdf, file=file) 
 # colors <- colorRampPalette(c('green', 'red'))(256)
 # fig<-levelplot(class$Sx,col.regions=colors, cuts=100,scales=list(
 #   x=list(rot=90)
 # ))
 # print(fig)
 # graphics.off()
  return(class)
}

calclus.spec <- function(data, index, kl , Sx,tags,   diss=FALSE, debug=FALSE, FUNCluster, FUNTest,...)
  #1=pknng, 2 = rbf, 3= local scaling,
{
  
  
  tags[index]<-kl
  t1 <- proc.time()

  sc<-FUNCluster(data=data,index=index,diss=diss,debug, ...)
  print("clustering")
  print(proc.time() - t1)
  if(!sc$error){
    
    left <- index[sc$clusters==1]  
    right<- index[sc$clusters==2]

    ind<-sort(c(left,right))   
    print(length(left))
#     plot(data[ind,],col=sc$clusters)
 
    if(kl==1) {
      
     Sx<-sc$Sx
     
    }
    

    tags[left] <- 2*kl 
    tags[right] <- 2*kl+1
    error<-0
    t1 <- proc.time()
    g<-FUNTest(data[ind,],sc,...)
    print("validation")
    print(proc.time() - t1)
  
  }
  else #spectral error
  {
    print("Error")  
    error <-1
  }


 if(!error && g$cond ) 
  { print(length(left))
   if(kl>1) 
    {
      innd<-1:length(ind)
      il<-innd[sc$clusters==1]  
      ir<-innd[sc$clusters==2] 
      Sx[left,left]<-sc$Sx[il,il]
      Sx[right,right]<-sc$Sx[ir,ir]
#       colors <- colorRampPalette(c('green', 'red'))(256)
#       fig<-levelplot(Sx,col.regions=colors,scales=list(
#         x=list(rot=90)
#       ))
#       print(fig)
    }
    print("true")
    kll<-kl
    klr<-kl
    if (length(left)>(3))
    {
      bl <- calclus.spec (data,  index=left, 2*kll ,Sx=Sx ,tags,  diss=diss,debug=debug,FUNCluster,FUNTest,...) 
    #  final <-bl
      tags[left]=(bl[[2]])[left]
      Sx[left,left]=(bl[[1]])[left,left]
      rm(bl)
    }else
    { 
      tags[left]<-kl      
    }
    if (length(right)>(3))
    {      
      br <- calclus.spec (data,  index=right,  (2*klr+1) ,Sx=Sx,tags,diss=diss, debug=debug,FUNCluster,FUNTest,...)
      tags[right]=(br[[2]])[right]
      Sx[right,right]=(br[[1]])[right,right]
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

  final<-list(Sx,tags)
  names(final)<-list("Sx", "tags")
  return(final)
}









