\name{test_gap}
\alias{test_gap}

\title{A test for Cluster tendency based on Gap Statistic}
\description{
test_gap() calculates a goodness of clustering measure, the gap statistic for 1 and 2 clusters.
For each number of clusters k=1:2,
it compares log(W(k)) with E*[log(W(k))]
where the latter is defined via bootstrapping, i.e. 
simulating from a reference distribution.

}

\usage{
test_gap(data, labels, NumRef=100, sigmas=1,
k=7, sigma=1,  method=2, debug=FALSE, metric=2)
}

\arguments{
\item{data}{matrix of data.}
\item{labels}{a vector of length n = nrow(data) of integers in 1:2 determining the clustering of the data.}
\item{NumRef}{integer, number of Monte Carlo samples.}
\item{debug}{Logical, determining if "progress" output should be printed. }
\item{method}{Integer indicating which method will be used to compute the 
similarity matrix: 1 for RBF with pknng, 2 for RBF with euclidean and 3
for Local Scaling }
\item{sigma}{.}
\item{sigmas}{.}
\item{k}{.}
\item{metric}{.}
}
\value{
a list with the next components
\item{Gap_info}{a matrix with 2 rows and 4 columns, named "OlogW", "E.logW", "Gap", and "Gap.sd", 
where Gap = E.logW - ologW, and Gap.sd corresponds to the standard error of gap, Gap.sd[k]=s[k], 
where s[k] := sqrt(1 + 1/B) sd^*(Gap[]), and sd^*() is the standard deviation of the simulated ("bootstrapped") gap values.}
\item{cond}{Boolean value result of the evaluation of the condition Gap[1]<Gap[2]-Gap.sd[2].}
\item{gap}{Distance between Gap[1] and Gap[2]-Gap.sd[2].}

}

\seealso{
\code{\link{pknng}}.
 }
\examples{
## 

}


\keyword{custering}
