\name{boot.vc}
\alias{boot.vc}
\title{Bootstrap confidence intervals for variance components}
\description{Provides a bootstrap confidence interval (over loci) for sums of the different variance components (equivalent to gene diversity estimates at the different levels), and the  derived F-statistics, as suggested by Weir and Cockerham (1984). Will not run with less than 5 loci. Raymond and Rousset (199X) points out shortcomings of this method.}
\usage{
boot.vc(levels=levels,loci=loci,diploid=TRUE,nboot=1000,quant=c(0.025,0.5,0.975))
}
\arguments{
\item{levels}{a data frame containing the different levels (factors) from the outermost (e.g. region) to the innermost before the individual}
\item{loci}{a data frame containing the different loci}
\item{diploid}{Specify whether the data are coming from diploid or haploid organisms (diploid is the default)}
\item{nboot}{Specify the number of bootstrap to carry out. Default is 1000}
\item{quant}{Specify which quantile to produce. Default is c(0.025,0.5,0.975) giving the percentile 95\% CI and the median}
}
\value{
\item{boot}{a data frame with the bootstrapped variance components. Could be used for obtaining bootstrap ci of statistics not listed here.}
\item{res}{a data frame with the bootstrap derived statistics. H stands for gene diversity, F for F-statistics}
\item{ci}{Confidence interval for each statistic.}
 }
\references{
Raymond M and Rousset F, 1995. An exact test for population differentiation. Evolution. 49:1280-1283

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Weir BS and Cockerham CC, 1984. Estimating F-statistics for the analysis of population structure. Evolution 38:1358-1370.

}

\seealso{\code{\link{varcomp.glob}}.}
\examples{
#load data set
data(gtrunchier)
boot.vc(gtrunchier[,c(1:2)],gtrunchier[,-c(1:2)],nboot=100)
}
\keyword{univar}
