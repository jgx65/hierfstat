\name{pop.freq}
\alias{pop.freq}
\title{Allelic frequencies}
\description{Estimates allelic frequencies for each population and locus}
\usage{pop.freq(dat,diploid=TRUE)}
\arguments{
\item{dat}{a data frame where the first column contains the population to which the different individuals belong, and the following columns contain the genotype of the individuals -one locus per column- }
\item{diploid}{specify whether the data set consists of diploid (default) or haploid data}
}

\value{
\item{}{A list containing allele frequencies. Each element of the list is one locus.
For each locus, Populations are in columns and alleles in rows}
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
pop.freq(gtrunchier[,-2])
}
\keyword{univar}
