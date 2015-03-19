\name{genet.dist}
\alias{genet.dist}
\title{Classical genetic distances estimation}
\description{
Estimates one of several genetic distances among all pairs of populations.
}
\usage{
genet.dist(dat,diploid=TRUE,method=\dQuote{Dch})
}
\arguments{
  \item{dat}{A data frame containing population of origin as the first column and multi-locus genotypes in following columns}
  \item{diploid}{whether the data is from a diploid (default) or haploid organism.}
  \item{method}{One of \dQuote{Dch},\dQuote{Da},\dQuote{Ds},\dQuote{Fst},\dQuote{Dm},\dQuote{Dr},\dQuote{Cp} or \dQuote{X2}, all described in Takezaki and Nei (1996)}
}

\details{the method argument specify which genetic distance to use, among eight, all briefly described in Takezaki and Nei (1996)

\item{\dQuote{Dch}}{By default, Cavalli-Sforza and Edwards Chord distance (eqn 6 in the reference) is returned. 
This distance is used as default as Takezaki & Nei (1996) found that it was the best to retrieve the relation among samples}
\item {\dQuote{Da}}{This is Nei's et al genetic distance (eqn 7), performing nearly as well as \dQuote{Dch}}
\item{\dQuote{Ds}}{Nei's standard genetic distance (eqn 1). Increases linearly with diverence time but has larger variance }
\item{\dQuote{Fst}}{Latter's and also approximately Reynolds et al Genetic distance (eqn 3)}
\item{\dQuote{Dm}}{Nei's minimum distance (eqn 2)}
\item{\dQuote{Dr}}{Rogers's distance (eqn 4)}
\item{\dQuote{Cp}}{Prevosti et al's distance (eqn 5)}
\item{\dQuote{X2}}{Sanghvi's distance (eqn 8)}
}
\value{
A matrix of pairwise genetic distance
}
\references{takezaki & Nei (1996) Genetic distances and reconstruction of Phylogenetic trees from microsatellite DNA. Genetics 144:389-399. \url{http://www.genetics.org/content/144/1/389}
}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

\seealso{\code{\link{fst.pp}}}

\examples{
data(gtrunchier)
genet.dist(gtrunchier[,-1])
genet.dist(gtrunchier[,-1],method="Dr")
}

\keyword{distance}
