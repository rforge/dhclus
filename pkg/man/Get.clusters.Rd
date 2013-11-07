\name{Get.clusters}
\alias{Get.clusters}
\title{Spectral clustering with 2 centers}
\description{
Get.clusters() applies spectral clustering with 2 center

}

\usage{
Get.clusters<-function(data, index, diss=FALSE, debug=FALSE,method=2)
}

\arguments{
\item{data}{matrix of data}
\item{index}{integer vector indicating which element of data will be used}
\item{diss}{Logical, Indicate if the data matrix has data information or disimilarity information.}
\item{debug}{Logical, determining if "progress" output should be printed. }
\item{method}{Integer indicating which method will be used to compute the 
similarity matrix: 1 for RBF with pknng, 2 for RBF with euclidean and 3
for Local Scaling
}

}
\seealso{
\code{\link{pknng}}.
 }
\examples{
## 

}


\keyword{custering}
