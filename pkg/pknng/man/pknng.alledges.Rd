\name{pknng.alledges}

\alias{pknng.alledges}

\title{Connects a set of point using the minimum number of edges}

\description{
  Given a dataset this method calculates the knn-graph. The method returns a list that represents a graph.
  This method connects each point in every group of neighbours to all the rest of the points not belonging 
  to the current group using a penalized metric.
}

\usage{
pknng.alledges (X,k,method="euclidean")
}

\arguments{
  \item{X}{Dataset. Cols representes dimensions and rows the number of points.}
  \item{k}{k is the number of neighbours.}
  \item{method}{the metric/similarity used, currently available: "euclidean,
  "mi" mutual information, "cor" correlation coefficient.}
}

\value{
  \code{pknng.alledges} returns a list containing
  \item{G}{Matrix with minimum distances between points. (The method uses 
  Floyd algorithm to obtain all pairs of distances).}
  \item{mu}{Mean value of k-nearest neighbours distance.}
  \item{psi}{After connecting the groups of knn, the new Mean distance value.}

}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.minspan}} and \code{\link{pknng.allsubgraphs}}.}

\examples{
## Calculates knn-graph (k=3) using euclidean metric, penalization and nearest point for groups union

library(pknng)
library(cluster)
rm(list=ls())

data(croissants,package="pknng")

# number of clusters
centers<- length(unique(croissants2d[,3]))

plot(croissants2d[,1],croissants2d[,2],col=croissants2d[,3],cex=0.25)
iso.croiss<-pknng.alledges (croissants2d[,-3],3,method="euclidean")

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
