\name{allele.count}
\alias{allele.count}
\title{Allelic counts}
\description{Counts the number of copies of the different alleles at each locus and population}
\usage{allele.count(data,diploid=TRUE)}
\arguments{
\item{data}{A data frame containing the population of origin in the first column and the genotypes in the following ones}
\item{diploid}{Whether the data are from diploid individuals}
}

\value{
A list of tables, --each with np (number of populations) columns and nl (number of loci) rows-- of the count of each allele
}

%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
allele.count(gtrunchier[,-2])
}
\keyword{univar}
