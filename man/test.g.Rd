\name{test.g}
\alias{test.g}
\title{Tests the significance of the effect of level on genetic differentiation}
\description{Tests the significance of the effect of \kbd{level} on genetic differentiation}
\usage{
test.g(data = data, level, nperm = 100,...)
}
\arguments{
\item{data}{a data frame containing the genotypes for the different loci}
\item{level}{A vector containing the assignment of each observation to its level}
\item{nperm}{The number of permutations to carry out for the test}
\item{...}{Mainly here to allow passing \kbd{diploid=FALSE} if necessary}
}
\value{
\item{g.star}{A vector containing all the generated g-statistics, the last one beeing the observed}
\item{p.val}{The p-value associated with the test}
}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\examples{
data(gtrunchier)
attach(gtrunchier)
test.g(gtrunchier[,-c(1,2)],Locality)
}
\keyword{nonparametric}
