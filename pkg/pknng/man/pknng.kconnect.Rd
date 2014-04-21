\name{pknng.kconnect}

\alias{pknng.kconnect}

\title{Connects a set of point using k neigbours }

\description{
  Given a dataset this method calculates the knn-graph. The method returns a list that represents a graph.
}

\usage{
pknng.kconnect (X,k,method="euclidean",conntype="minspan",union="near",penalize=1)
}

\arguments{
  \item{X}{Matrix. Cols representes dimensions and rows the number of points.}
  \item{k}{Integer. k is the number of neighbours.}
  \item{method}{Type of metric/similarity used, currently available: "euclidean,
  "mi" mutual information, "cor" correlation coefficient.}
  \item{conntype}{Type of connection algorithm. Currently "minspan" and "allsubgraphs". See \code{\link{pknng.allsubgraphs}}}
  \item{union}{Subtype of "allsubgraphs" connection algorithm. Currently "union" and "medoid". See \code{\link{pknng.allsubgraphs}}}
  \item{penalize}{Type. penalize = 1 use penalization, penalize = 0 does not use penalization.}
}

\value{
  \code{pknng.kconnect} returns a list containing
  \item{i}{Integer vector. Source point which means that an edge goes from i to j.}
  \item{k}Integer vector. Destination point which means that an edge goes from i to j.}
  \item{s}{Real vector. Represents distance from point i to point j.}

}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.getmink}} and \code{\link{pknng.allsubgraphs}} and \code{\link{pknng.apcluster}}.}

\examples{
## Calculates knn-graph using euclidean metric, penalization and return a sparse matrix
pknng.kconnect (X)

## Calculates knn-graph using mutual information, penalization and return a sparse matrix
pknng.kconnect (X,method="mi")
}
\keyword{graphs} 
