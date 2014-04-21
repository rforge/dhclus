\name{pknng.report}

\alias{pknng.report}

\title{Connects a set of point using k neigbours }

\description{
  Given a dataset this method calculates the knn-graph. The method returns a list that represents a graph.
}

\usage{
pknng.report (class,pclass)
}

\arguments{
  \item{class}{Integer vector. True data class.}
  \item{pclass}{Integer vector. Predicted class.}
}

\value{
  \code{pknng.report} returns a list containing
  \item{ac}{Real value. Accuracy}
  \item{vector.pclass}{Real vector representing the accuracy for each predicted class refered to the true class.}
  \item{ctable}{Matrix. Confusion matrix}

}

\author{Ariel Baya \email{abaya@ifir.edu.ar}, Pablo Granitto.}

\seealso{See also  \code{\link{pknng.getmink}}.}

\examples{
## Calculates accuracy of predicted vector pclass_a with respect to class_a
pknng.report (class_a,pclass_a)

}
\keyword{graphs}  
