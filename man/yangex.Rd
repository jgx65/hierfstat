\name{yangex}
\docType{data}
\alias{yangex}
\title{Example data set from Yang (1998) appendix}
\description{Reproduce the example data set used in Yang's paper appendix. 
The genotype (column genot) is invented}
\usage{
data(exhier)
}
\value{
\item{pop}{outermost level}
\item{spop}{sub pop level}
\item{sspop}{sub sub pop level}
\item{genot}{dummy diploid genotype}

}
\references{
Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52(4):950-956
}

\examples{
data(yangex)
varcomp(yangex)
#the k matrix should be the same as matrix (A2) in Yang's appendix, p. 956
}
\keyword{datasets}
