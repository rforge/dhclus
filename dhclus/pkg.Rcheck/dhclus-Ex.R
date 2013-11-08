pkgname <- "dhclus"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('dhclus')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Get.clusters")
### * Get.clusters

flush(stderr()); flush(stdout())

### Name: Get.clusters
### Title: Spectral clustering with 2 centers
### Keywords: custering

### ** Examples

## 




cleanEx()
nameEx("dhclus")
### * dhclus

flush(stderr()); flush(stdout())

### Name: dhclus
### Title: Divisive hierarchical clustering algorithm with automatic
###   determination of components.
### Keywords: custering

### ** Examples

## Using default spetral clustering


data<-cbind(rnorm(50,0,1),rnorm(100,0,1),rep(1,100))
data<-rbind(data, cbind(rnorm(50,0,1)+6.5,rnorm(100,0,1),rep(2,100)))
data<-rbind(data, cbind(rnorm(20,0,1)-5,rnorm(100,0,1)-5,rep(3,100)))

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




cleanEx()
nameEx("test_gap")
### * test_gap

flush(stderr()); flush(stdout())

### Name: test_gap
### Title: A test for Cluster tendency based on Gap Statistic
### Keywords: custering

### ** Examples

## 




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
