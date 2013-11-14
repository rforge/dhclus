\name{mst2Path.Diss}
\alias{mst2Path.Diss}

\title{MST based distance matrix}
\description{


}

\usage{
mst2Path.Diss(X)
}

\arguments{
\item{X}{matrix of distances of data .}

}
\value{
a matrix of MST based distances

}


\examples{
## 
library(icav)
data(fractal)
d<-mst2Path.Diss(as.matrix(dist(datos)))
}



