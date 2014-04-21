\name{pknng.cormetric}

\alias{pknng.cormetric}

\title{Calculates correlation coefficient (fast)}

\description{
  Calculates correlation coefficient.
}

\usage{
pknng.cormetric (data))
}

\arguments{
  \item{data}{Matrix. Cols represents the number of points and rows the dimensions of the data.}}

\value{
  \code{pknng.cormetric} Matrix of correlation coefficients
}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.eucmetric}} and \code{\link{pknng.mutualmetric}}.}

\examples{
## Calculates correlation coefficients matrix.
pknng.cormetric (Data)
}
\keyword{graphs} 
