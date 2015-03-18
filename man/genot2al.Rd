\name{genot2al}
\alias{genot2al}
\title{Separates diploid genotypes in its constituant alleles}
\description{
Separates the input vector of diploid genotypes in two vectors each containing one allele, and returns a vector of length \kbd{2*length(y)} with the second part being the second allele
}
\usage{
genot2al(y)
}
\arguments{
\item{y}{the diploid genotypes at one locus}
}
\value{
returns a vector of length \kbd{2*length(y)}, with the second half of the vector containing the second alleles}
\references{
Goudet J. (2004). A library for R to compute and test variance components and F-statistics. In Prep
}
\author{Jerome Goudet,
DEE, UNIL, 
CH-1015
Lausanne
Switzerland 

\email{jerome.goudet@unil.ch}
}
\seealso{\code{\link{varcomp}}.}
\examples{
data(gtrunchier)
genot2al(gtrunchier[,4])
}
\keyword{manip}
