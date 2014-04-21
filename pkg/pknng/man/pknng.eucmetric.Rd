\name{pknng.eucmetric}

\alias{pknng.eucmetric}

\title{Calculates euclidean distance matrix}

\description{
  Calculates euclidean distance matrix faster than dist(method="euclidean"). 
}

\usage{
pknng.eucmetric (data))
}

\arguments{
  \item{data}{Matrix. Cols represents the number of points and rows the dimensions of the data.}}

\value{
  \code{pknng.eucmetric} Euclidean distance matrix
}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.cormetric}} and \code{\link{pknng.mutualmetric}}.}

\examples{
## Calculates euclidean distance matrix.
pknng.eucmetric (Data)
}
\keyword{graphs} 
 
