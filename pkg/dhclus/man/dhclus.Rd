\name{dhclus}
\alias{dhclus}
\title{Divisive hierarchical clustering algorithm with automatic determination of components.}
\description{
Generic function for clustering data using hierarchical strategies. 
By default uses spectral clustering with RBF kernel 
as partitional clustering strategy in each level. 

}

\usage{
dhclus(data, diss=FALSE,  FUNCluster=Get.clusters, FUNTest=test, ...)
}

\arguments{
\item{data}{matrix of data}
\item{diss}{Logical, indicates if the data matrix has data information or disimilarity information.}
\item{FUNCluster}{
a function which accepts as first argument a (data) matrix, second argument an index vector
indicating which element of data will be used, and returns a list with a component named   clusters which is a vector of length n = length(index) of integers in 1:2 determining the clustering or grouping of the n observations and a component named error which is a boolean
value indicates whether an error has occurred to group items.
}
\item{FUNTest}{a function which accepts as first argument a (data) matrix, second
argument a vector of length n = nrow(data) of integers in 1:2 determining the clustering of the data, and returns a list with a component named cond which is a boolean value 
determines that the given clustering must be accepted or rejected.}
\item{...}{ optionally further arguments for FUNcluster() and FUNTest().}
}
\value{
a list whith clustering information
\item{tabs}{}  
\item{ltabs}{} 
\item{labelsTree}{} 
\item{as}{} 
\item{ssb}{} 
\item{tags_sil}{} 
\item{tags_bss}{} 
\item{heigth}{} 
}

\seealso{
\code{\link{Get.clusters},\link{test_gap}}.
 }
\examples{
## Using default spetral clustering


data<-cbind(rnorm(70,0,1),rnorm(70,0,1),rep(1,70))
data<-rbind(data, cbind(rnorm(50,0,1)+6.5,rnorm(50,0,1),rep(2,50)))
data<-rbind(data, cbind(rnorm(30,0,1)-5,rnorm(30,0,1)-5,rep(3,30)))

clusters<- dhclus(data,debug=FALSE)
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
MyTest<-function(data,result,...)
{
  return(test_gap(data,labels=result$clusters,NumRef=50,metric=3,...))
}


clusters<- dhclus(data, FUNCluster=Myclustering, FUNTest=MyTest) 

}
\keyword{custering}
