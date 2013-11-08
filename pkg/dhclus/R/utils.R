

#Distsum: Computes the sum of pairwise distances

Distsum<-function(dist){
  d<-as.matrix(dist)
  D<-sum(rowSums(d))
  return(0.5*D)
}

#Distsum: Computes the sum of squared pairwise distances
Distsum2<-function(dist){
  d<-as.matrix(dist)
  D<-sum(rowSums(d^2))
  return(0.5*D)
}



# The pooled within-cluster sum of squares around the cluster means
# If diss is TRUE data are treated as a distance matrix

withinsum2<-function(data,labels, diss=FALSE){
  data<-as.matrix(data)
  
  if(length(labels) == dim(data)[1])
  {
    index<-1: dim(data)[1]
    labels<-relabel(labels) #in case labels not start with 1,2,etc
    l<-length(unique(labels))
    w<-vector(length=l)
    w<-sapply(1:l, function(i) {
      ind<-index[labels==i]
      n<-length(ind)
      if(diss==FALSE) 
        D<-0.5*(sum(dist(data[ind,])^2))
      else  
        D <- Distsum2(data[ind,ind])
      return(D/n)
    })
  } else  stop("length Labels must be the same of data ")
  return(w)
}


bss<-function(data,labels, diss=FALSE){
  data<-as.matrix(data)
  k<-length(unique(labels))
  w1<-withinsum2(data,rep(1,length(labels)),diss=diss)
  wk<-sum(withinsum2(data,labels,diss=diss))
  w<-w1-wk
  return(as.double((k*wk))/(w))
}

##functions for trees

max_level<-function(class){
  uni<-unique(class)
  return(max(floor(log2(uni))))
}

possibles<-function(i)
{
  return(2^i )
}

childrens<-function(i,max,j)
{
  if(i>max) out<-NULL
  if(i==max) out<-c(2*(2^i+j),2*(2^i+j)+1)
  if(i<max) out<-c(2*(2^i+j),2*(2^i+j)+1,childrens( i=2*(2^i+j), max,j),childrens( i=(2*(2^i+j)+1), max,j))

  return(out)
}

labels_level<-function(class,i){
  ##must be i>1
  index<-1:length(class)
  labs<-class
  c<-possibles(i) 
  
  ls<-((1:c)+2^i)-1

  m<-matrix(0,nrow=length(class),ncol=c/2)
  j<-1
  while(length(ls)>1){
    m[,j]<-labs
    lmax<-max(ls)
    cl<-c(lmax,lmax-1)
    ls<-ls[!is.element(ls,cl)]

    ind<-index[is.element(labs,cl)]    
    if(length(ind)>0){
      m[ind,j]<-floor((labs[ind])/2)
      labs<-m[,j]  
      j<-j+1}
  }
  ou<-cbind(m[,(j-1):1])
  
  return(ou)
}

##############3
first_nonzero<-function(v,k){
  flag<-0
  while(v[k]==0){
    k<-k+1
    if(k>=length(v)){ flag<-1;break;}
  }
  if(flag==1) return (1)
  else return (v[k])
}



relabel<-function(etiquetas)
{
  cl1<-rep(0,length(etiquetas))
  uni1<-unique(etiquetas); 
  for(i in 1:length(uni1)) {cl1[etiquetas==uni1[i]]<-i}
  return(cl1)
}



