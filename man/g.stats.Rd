\name{g.stats}
\alias{g.stats}
\title{Calculates likelihood-ratio G-statistic on contingency table}
\description{Calculates the likelihood ratio G-statistic on a contingency table of alleles at one locus X sampling unit. The sampling unit could be any hierarchical level}
\usage{
g.stats(data,diploid=TRUE)}
\arguments{
\item{data}{a two-column data frame. The first column contains the sampling unit, the second the genotypes}
\item{diploid}{Whether the data are from diploid (default) organisms}
}
\value{
\item{obs}{Observed contingency table}
\item{exp}{Expected number of allelic observations}
\item{X.squared}{The chi-squared statistics, \eqn{\sum\frac{(O-E)^2}{E}}{Sum (O-E)^2/E}}
\item{g.stats}{The likelihood ratio statistics, \eqn{2 \sum(O \log(\frac{O}{E}))}{2 Sum(O Log(O/E))}}
 }
\references{
Goudet J., Raymond, M., DeMeeus, T. and Rousset F. (1996) Testing differentiation in diploid populations. Genetics. 144: 1933-1940

Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

Petit E., Balloux F. and Goudet J.(2001) Sex-biased dispersal in a migratory bat: A characterization using sex-specific demographic parameters. Evolution 55: 635-640. 

}
\author{Jerome Goudet,
DEE, UNIL, 
CH-1015
Lausanne
Switzerland 

\email{jerome.goudet@unil.ch}
}
\seealso{\code{\link{g.stats.glob}}.}
\examples{
data(gtrunchier)
attach(gtrunchier)
g.stats(data.frame(Patch,L21.V))
}
\keyword{univar}
