\name{pknng.apcluster}

\alias{pknng.apcluster}

\title{Affinity propagation clustering method}

\description{
	Affinity propagation clustering method.
}

\usage{
pknng.apclustert (data,i,j,datalen,n,preference,damp,maxiter,conviter=50,noise=0)
}

\arguments{
  \item{data}{Real vector. Represents distance from point i to point j.}
  \item{i}{Integer vector. Source point which means that an edge goes from i to j.}
  \item{j}{Integer vector. Destination point which means that an edge goes from i to j.}
  \item{datalen}{Integer. Length of data, i, j. The three of them have to agree. This represents the number of edges in the graph.}
  \item{n}{Integer. Number of point in the the graph.}
  \item{preference}{Real. Preference for affinity clustering algorithm. Read the paper for more information.}
  \item{damp}{Real. Damping factor (0 < damp <= 1). Read the paper for more information.}
  \item{maxiter}{Integer. Maimum number of iteration.}
  \item{conviter}{Integer. Minimum number of iteration to evaluate convergence.}
  \item{noise}{Type. If noise=1 affinity adds some minimum noise to data (see paper for datails), 
   if noise=0 no noise added. If noise = 1 data IS CHANGED to data+noise.}
}

\value{
  \code{pknng.apcluster} returns a list containing
  \item{index}{Integer vector. The index of this vector refers to the point id and the actual 
  value in that position to its center. So this vector contains the centers id corresponding to 
  the index of this vector.}
  \item{param}{Vector. param[1] = Number of identified clusters, param[2] = Fitness (net similarity), 
  param[3] = Similarities of data points to exemplars, param[4] = Preferences of selected exemplars, 
  param[5] = Number of iterations}
  \item{class}{Vector. The classes assigned to the data.}

}

\details{
  C-implementation of the affinity propagation clustering algorithm. See  BJ Frey and D Dueck,
  Science 315, 972-976, Feb 16, 2007, for a full description of the algorithm. Copyright 2007, 
  BJ Frey and Delbert Dueck. This software may be freely used and distributed for non-commercial 
  purposes. Based on the code available at http://www.psi.toronto.edu/affinitypropagation/
  modified by Ariel E. Baya.
}

\references{ 
Frey B.J.,Dueck D.\cr
\emph{Clustering by passing messagesbetween data point.}\cr
Science, vol. 315, pp. 972-976, Feb 16 2007\cr
\url{http://www.psi.toronto.edu/affinitypropagation/}
}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.kconnect}}.}

\examples{
## Calculates knn-graph (k=3) using euclidean metric, penalization and nearest point for groups union
pknng.DistanceMatrix.one.totherest (X,3)

}
\keyword{graphs}  
