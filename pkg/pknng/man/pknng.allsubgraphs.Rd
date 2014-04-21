\name{pknng.allsubgraphs}

\alias{pknng.allsubgraphs}

\title{Connects a set of point using the minimum number of edges}

\description{
  Given a dataset this method calculates the knn-graph. The method returns a list that represents a graph.
}

\usage{
pknng.allsubgraphs (X,k,method="euclidean",union="near",penalize=1)
}

\arguments{
  \item{X}{Matrix. Cols representes dimensions and rows the number of points.}
  \item{k}{Integer. k is the number of neighbours.}
  \item{method}{Type of metric/similarity used, currently available: "euclidean,
  "mi" mutual information, "cor" correlation coefficient.}
  \item{union}{Type that specifies how to connetc the k-nn,"near" connects from one group of neighbours
  to all the rest throgh the nearest points of each two groups, while "medoid" connects through the 
  medoids (centroid) of each two groups.}
  \item{penalize}{Type. penalize = 1 use penalization, penalize = 0 does not use penalization.}
}

\value{
  \code{pknng.allsubgraphs} returns a list containing
  \item{G}{Matrix of minimum distances between points. (The method uses 
  Floyd algorithm to obtain all pairs of distances).}
  \item{mu}{Real value. Mean value of k-nearest neighbours distance.}
  \item{psi}{Real value. After connecting the groups of knn, the new Mean distance value.}

}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.minspan}} and \code{\link{pknng.alledges}}.}

\examples{
## Calculates knn-graph (k=3) using euclidean metric, penalization and nearest point for groups union

library(pknng)
library(cluster)
rm(list=ls())

data(croissants,package="pknng")

# number of clusters
centers<- length(unique(croissants2d[,3]))

plot(croissants2d[,1],croissants2d[,2],col=croissants2d[,3],cex=0.25)
iso.croiss<-pknng.allsubgraphs(croissants2d[,-3],3,method="euclidean",union="near",penalize=1)

dend<-hclust(as.dist(iso.croiss$G),method="single")
h<-1:20
dend.clust<-dend.error<-0.0
for(ii in 1:length(h)){
dend.l<-cutree(dend,k=h[ii])
dend.clust[ii]<-nlevels(factor(dend.l))
dend.error[ii]<-pknng.report(croissants2d[,3],dend.l)$ac
}

error.dend <-max(dend.error)
dend.clusters<-which.max(dend.error)
rm(dend.l);

pam.c<-pam(as.dist(iso.croiss$G),k=centers,diss=T)
error.pam<-pknng.report(croissants2d[,3],as.vector(pam.c$clustering))$ac
rm(pam.c)

print(paste("Clusters recognized by slink dendogram: ",dend.clusters," accuracy: ",error.dend))
print(paste("Accuracy of PAM with ",centers," centers : ",error.dend))

}
\keyword{graphs} 
