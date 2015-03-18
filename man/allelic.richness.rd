\name{allelic.richness}
\alias{allelic.richness}
\title{Estimates allelic richness}
\description{Estimates allelic richness, the rarefied allelic counts, per locus and population}
\usage{allelic.richness(data,min.n=NULL,diploid=TRUE)}
\arguments{
\item{data}{A data frame, with as many rows as individuals. The first column contains the population to which the individual belongs, the following to the different loci}
\item{min.n}{The number of alleles down to which the number of alleles should be rarefied. The default is the minimum number of individuals genotyped (times 2 for diploids)}
\item{diploid}{a boolean specifying wether individuals are diploid (default) or haploid}
}

\value{
\item{min.all}{The number of alleles used for rarefaction}
\item{Ar}{A table with as many rows as loci and columns as populations containing the rarefied allele counts}
} 
\references{
El Mousadik A. and  Petit R.J. (1996) High level of genetic differentiation for allelic richness among populations of the argan tree argania spinosa skeels endemic to Morocco. Theoretical and Applied Genetics, 92:832-839

Hurlbert S.H. (1971) The nonconcept of species diversity: a critique and alternative parameters. Ecology, 52:577-586

Petit R.J.,  El Mousadik A. and Pons O. (1998) Identifying populations for conservation on the basis of genetic markers. Conservation Biology, 12:844-855
}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

\examples{
data(gtrunchier)
allelic.richness(gtrunchier[,-1])
}
\keyword{univar}
