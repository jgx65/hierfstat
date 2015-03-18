\name{test.within}
\alias{test.within}
\title{Tests the significance of the effect of inner.level on genetic differentiation within blocks defined by outer.level}
\description{Tests the significance of the effect of inner.level on genetic differentiation within blocks defined by outer.level}
\usage{
test.within(data, within, test.lev, nperm, ...)
}
\arguments{
\item{data}{a data frame containing the genotypes for the different loci}
\item{within}{A vector containing the units in which to keep the observations}
\item{test.lev}{A vector containing the units from which to construct the contingency tables}
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
#tests whether the patch level has a significant effect on genetic structure
test.within(gtrunchier[,-c(1,2)],within=Locality,test.lev=Patch)
}
\keyword{nonparametric}
