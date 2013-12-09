\name{Get.clusters}
\alias{Get.clusters}
\title{Spectral clustering with 2 centers}
\description{
Get.clusters() applies spectral clustering with 2 center

}

\usage{
Get.clusters(data, index, diss=FALSE, debug=FALSE, centers=2, method=2,...)
}

\arguments{
\item{data}{matrix of data}
\item{index}{integer vector indicating which element of data will be used}
\item{diss}{Logical, Indicate if the data matrix has data information or disimilarity information.}
\item{debug}{Logical, determining if "progress" output should be printed. }
\item{centers}{integer, the number of clusters. }
\item{method}{Integer indicating which method will be used to compute the 
similarity matrix: 1 for RBF with pknng, 2 for RBF with euclidean and 3
for Local Scaling.
}
\item{...}{ optionally further arguments.}
}
\value{
a list with the next components
\item{clusters}{ vector of length n = length(index) of integers in 1:2 determining the clustering of the data.}
\item{val}{.}
\item{sigma}{value}
\item{sigmas}{.}
\item{error}{boolean value indicates if an error whether an error has occurred to group items}

}

\seealso{
\code{\link{pknng}}.
 }
\examples{
## 

}


\keyword{custering}
