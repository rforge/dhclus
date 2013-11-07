\name{dhclus}
\alias{dhclus}
\title{Divisive hierarchical clustering algorithm with automatic determination of components.}
\description{


}

\usage{
dhclus<- function(data, diss=FALSE,debug=FALSE,FUNCluster=Get.clusters, FUNTest=test,...)
}

\arguments{
\item{data}{matrix of data}
\item{diss}{Logical, Indicate if the data matrix has data information or disimilarity information.}
\item{debug}{Logical, determining if “progress” output should be printed. }
\item{FUNCluster}{}
\item{FUNTest}{}
}
\seealso{
\code{\link{Get.clusters}}.
 }
\examples{
## Using default spetral clustering


data<-cbind(rnorm(50,0,1),rnorm(100,0,1),rep(1,100))
data<-rbind(data, cbind(rnorm(50,0,1)+6.5,rnorm(100,0,1),rep(2,100)))
data<-rbind(data, cbind(rnorm(20,0,1)-5,rnorm(100,0,1)-5,rep(3,100)))

clusters<- dhclus(data)
plot(data, col=clusters$tags)

## Using spectral clustering with Local Scaling
clusters<- dhclus(data,debug=FALSE, method=3)

## Using spectral clustering with pknng and test-gap-pknng
clusters<- dhclus(data, debug=FALSE, method=1,metric=1)

##Using own created functions
# 

Myclustering<-function(data,index,...)
{
  res<-kmeans(data[index,],2)
  res$error<-FALSE
  res$clusters<-res$cluster
  return(res)
}

# Euclidean gap test
#MyTest<-function(data,result,...)
#{
#  return(test_gap(data,labels=result$cluster,NumRef=50,...))
#}


#clusters<- dhclus(data, FUNCluster=Myclustering, FUNTest=MyTest) 

}
\keyword{custering}
