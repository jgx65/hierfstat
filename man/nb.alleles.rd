\name{nb.alleles}
\alias{nb.alleles}
\title{Number of different alleles}
\description{Counts the number of different alleles at each locus and population}
\usage{nb.alleles(data,diploid=TRUE)}
\arguments{
\item{data}{A data frame containing the population of origin in the first column and the genotypes in the following ones}
\item{diploid}{whether individuals are diploid}
}

\value{
\item{}{A table, --with np (number of populations) columns and nl (number of loci) rows-- of the number of different alleles}
}
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
nb.alleles(gtrunchier[,-2])
}
\keyword{univar}
