\name{ind.count}
\alias{ind.count}
\title{individual counts}
\description{Counts the number of individual genotyped per locus and population}
\usage{ind.count(data)}
\arguments{
\item{data}{a data frame containing the population of origin in the first column and the genotypes in the following ones}
}

\value{
A table --with np (number of populations) columns and nl (number of loci) rows-- of genotype counts
}
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
ind.count(gtrunchier[,-2])
}
\keyword{univar}
