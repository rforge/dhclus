\name{pknng}

\alias{pknng}

\title{Connects a set of point using a K-nearest neighbours graph}

\description{
  Given a dataset this method calculates the k-nn graph returning a distance/similarity matrix. For a conectivity parameter the algorithm calculates all nearest neighbours graphs, 
then it connects them using a penalization. Finally it calculates all edges length using Dijkstra algorithm.
}

\usage{
pknng(X,k,diss=T,cte=3,mu="mean",method="euclidean",conn="one",penalize=1)
}

\arguments{
  \item{X}{Dataset where cols representes dimensions and rows the number of points or an n by n distance/similarity matrix.}

  \item{k}{k indicates the number of neighbours used to calculate all nearest neighbours graphs.}

  \item{diss}{Boolean value indicating if it is a distance/similarity matrix or not. Default value is TRUE.}

  \item{cte}{For penalize=2 or 3 it indicates the value of the lineal constant or the value of the polynomial exponent}

  \item{mu}{In the penalized form a penalized weigth (w) is obtained as w = w exp(w/mu) (exponencial penalization), w = wK (lineal penalization) and w = w (w/mu)^{K} (polynomic penalization). Values of mu can be: "mean", "median" , "1q" (first quartile), "3q" (thrid quartile)}

  \item{method}{if diss is FALSE. Method defines the metric/similarity used betweeen paterns, currently available: "euclidean", "mi.1" or "mi.2" mutual information, "corr" correlation coefficient.}
  
  \item{conn}{Method used for connection, "one" minimum spaning tree, "ttr" == connects each graph to all others using the closest points, "all" == conects every point from each graph to all other point not belonging to the same graph.}

  \item{penalize}{Integer type.  penalize=0 for non penalized conection, penalize=1 conects using exponencial penalization, penalize=2 conects using lineal penalization,penalize=3 conects 
using polynomic penalization.}
}

\value{
  \code{pknng} returns a distance/similarity matrix.
}

\author{Ariel Baya \email{baya@cifasis-conicet.gov.ar}, Pablo Granitto.}


\examples{
## Calculates knn-graph (k=3) using euclidean metric, penalization and MST for groups union
## See also test examples

library(pknng)
data(Three.Rings.LowNoise)

d.iso.a<-pknng(data[,-3],3,diss=FALSE,mu="mean",method="euclidean",conn="one",penalize=1)
HCmethod<-"average"

a<-hclust(as.dist(d.iso.a),method = HCmethod)
label.a<-cutree(a,k=5)
 split.screen(c(1,2))
 screen(1)
plot(data[,1],data[,2],col=data[,3],pch=16,cex=0.5)
title("True labels")

screen(2)
plot(data[,1],data[,2],col=label.a,pch=16,cex=0.5)
title("Clustering Result")
close.screen(all = TRUE) 
##see tests files in pknng/tests
}
\keyword{clustering}
