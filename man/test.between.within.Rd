\name{test.between.within}
\alias{test.between.within}
\title{Tests the significance of the effect of test.lev on genetic differentiation}
\description{Tests, using permutations of \kbd{rand.unit} within units defined by the vector \kbd{within} the significance of the contingency tables allele X (levels of \kbd{test.lev})}
\usage{
test.between.within(data, within, test.lev, rand.unit, nperm, ...) 
}
\arguments{
\item{data}{a data frame containing the genotypes for the different loci}
\item{within}{A vector containing the units in which to keep the observations}
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
data(yangex)
attach(yangex)
#tests for the effect of spop on genetic structure
test.between.within(data.frame(genot),within=pop,test=spop,rand=sspop)
}
\keyword{nonparametric}
