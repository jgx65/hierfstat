\name{genet.dist}
\alias{genet.dist}
\title{Classical genetic distances estimation}
\description{
Estimates one of several genetic distances among all pairs of populations.
}
\usage{
genet.dist(dat,diploid=TRUE,method="Dch")
}
\arguments{
  \item{dat}{A data frame containing population of origin as the first column and multi-locus genotypes in following columns}
  \item{diploid}{whether the data is from a diploid (default) or haploid organism.}
  \item{method}{One of \dQuote{Dch},\dQuote{Da},\dQuote{Ds},\dQuote{Fst},\dQuote{Dm},\dQuote{Dr},\dQuote{Cp} or \dQuote{X2}, all described in Takezaki and Nei (1996). Additionally \dQuote{Nei87} and \dQuote{WC84} return pairwise FSTs estimated following Nei (1987) \link{pairwise.neifst} and Weir & Cockerham (1984) \link{pp.fst} respectively}
}

\details{the method argument specify which genetic distance to use, among eight, all briefly described in Takezaki and Nei (1996)

\dQuote{Dch} By default, Cavalli-Sforza and Edwards Chord distance (eqn 6 in the reference) is returned. 
This distance is used as default since Takezaki & Nei (1996) found that it was the best to retrieve the relation among samples.

\dQuote{Da} This is Nei's et al genetic distance (eqn 7), performing nearly as well as \dQuote{Dch}

\dQuote{Ds} Nei's standard genetic distance (eqn 1). Increases linearly with diverence time but has larger variance 

\dQuote{Fst} Latter's and also approximately Reynolds et al Genetic distance (eqn 3)

\dQuote{Dm} Nei's minimum distance (eqn 2)

\dQuote{Dr} Rogers's distance (eqn 4)

\dQuote{Cp} Prevosti et al's distance (eqn 5)

\dQuote{X2} Sanghvi's distance (eqn 8) 

\dQuote{Nei87} see \link{pairwise.neifst}

\dQuote{WC84}  see \link{pairwise.WCfst}
}
\value{
A matrix of pairwise genetic distance
}
\references{Takezaki & Nei (1996) Genetic distances and reconstruction of Phylogenetic trees from microsatellite DNA. Genetics 144:389-399. \url{http://www.genetics.org/content/144/1/389}

Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press

Weir, B.S. and Cockerham C.C. (1984) Estimating F-Statistics for the Analysis of Population Structure 38:1358-1370 \url{http://www.jstor.org/stable/2408641}
}

\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

\seealso{\code{\link{pp.fst}} \code{\link{pairwise.neifst}}}

\examples{
data(gtrunchier)
genet.dist(gtrunchier[,-1])
genet.dist(gtrunchier[,-1],method="Dr")
}

\keyword{distance}
