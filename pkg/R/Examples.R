
#Clustering functions

Myclustering<-function(data,index,...)
{
  res<-kmeans(data[index,],2)
  res$error<-FALSE
  res$clusters<-res$cluster
  return(res)
}


###Test functions
MyTest<-function(data,result,...)
{
  return(test_gap(data,label=result$cluster,NumRef=50,...))
}