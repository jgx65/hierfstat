\name{test.between}
\alias{test.between}
\title{Tests the significance of the effect of test.lev on genetic differentiation}
\description{Tests the significance of the effect of test.lev on genetic differentiation}
\usage{
test.between(data, test.lev, rand.unit, nperm, ...)
}
\arguments{
\item{data}{a data frame containing the genotypes for the different loci}
\item{test.lev}{A vector containing the units from which to construct the contingency tables}
\item{rand.unit}{A vector containing the assignment of each observation to the units to be permutted}
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
#test whether the locality level has a significant effect on genetic structuring
test.between(gtrunchier[,-c(1,2)],test.lev=Locality,rand.unit=Patch)
}
\keyword{nonparametric}
