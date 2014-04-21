\name{pknng.getmink}

\alias{pknng.getmink}

\title{Connects a set of point using the minimum number of edges}

\description{
  Given a dataset this method calculates the knn-graph with the minimum K that results in a fully connected graph. 
The method returns a matrix representation of the graph.
}

\usage{
pknng.getmink (X,method="euclidean",return.matrix=FALSE,penalize=1))
}

\arguments{
  \item{X}{Matrix. Cols representes dimensions and rows the number of points.}
  \item{method}{Type of metric/similarity used, currently available: "euclidean,
 "mi" mutual information, "cor" correlation coefficient.}
  \item{return.matrix}{TRUE/FALSE. TRUE: after connecting the neighbours the method uses the floyd 
algorithm to obtaing all minimum distances between all pairs of points, FALSE: the method connects the k-neighbours and this groups of neighbours between them to obtain a fully connected graph leaving some distances undetermined (equal to -1).}
  \item{penalize}{Type. penalize = 1 use penalization, penalize = 0 does not use penalization.}
}

\value{
  \code{pknng.getmink} returns a list containing
  \item{G}{Matrix of distances between points}
  \item{k}{Integer value. Minimun k that fully connects this data}
}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.kconnect}}.}

\examples{
## Calculates knn-graph using euclidean metric, penalization and return a sparse matrix (does not use Floyd)
pknng.getmink (X)

## Calculates knn-graph using euclidean metric, penalization and return a matrix with distance to all pairs of points (uses Floyd)
pknng.getmink (X,return.matrix=TRUE)

## Calculates knn-graph using mutual information, penalization and return a sparse matrix (does not use Floyd)
pknng.getmink (X,method="mi")
}
\keyword{graphs}
