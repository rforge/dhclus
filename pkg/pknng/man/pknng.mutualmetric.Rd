\name{pknng.mutualmetric}

\alias{pknng.mutualmetric}

\title{Calculates mutual information }

\description{
  Calculates mutual information.
}

\usage{
pknng.mutualmetric (data,bins=0))
}

\arguments{
  \item{data}{Matrix. Cols represents the number of points and rows the dimensions of the data.}
  \item{bins}{Integer value. (bins*bins) is the number of bins used to construct the 2-d histogram. 
   If bins = 0, the number is calculated automatically.}}

\value{
  \code{pknng.mutualmetric} Matrix of correlation coefficients
}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.cormetric}} and \code{\link{pknng.eucmetric}}.}

\examples{
## Calculates mutual information matrix with automatic bins.
pknng.mutualmetric (Data)

## Calculates mutual information matrix with bins=49.
pknng.mutualmetric (Data,bins=7)
}
\keyword{graphs}  
