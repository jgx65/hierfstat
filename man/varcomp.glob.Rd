\name{varcomp.glob}
\alias{varcomp.glob}
\title{Estimate variance components and hierarchical F-statistics over all loci}
\description{Return multilocus estimators of variance components and F-statistics}
\usage{
varcomp.glob(levels=levels,loci=loci,diploid=TRUE)
}
\arguments{
\item{levels}{a data frame containing the different levels (factors) from the outermost (e.g. region) to the innermost before the individual}
\item{loci}{a data frame containing the different loci}
\item{diploid}{Specify whether the data are coming from diploid or haploid organisms (diploid is the default)}
}
\value{
\item{loc}{The variance components for each locus}
\item{overall}{The variance components summed over all loci}
\item{F}{a matrix of hierarchical F-statistics type-coefficients
 with the first line corresponding to \eqn{F_{(n-1)/n}},\eqn{F_{(n-2)/n}}...\eqn{F_{i/n}} and
 the diagonal corresponding to \eqn{F_{(n-1)/n}}, \eqn{F_{(n-2)/(n-1)}},\eqn{F_{i/2}}}
 }
\references{
Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52(4):950-956

Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186
}
\author{Jerome Goudet
DEE, UNIL, 
CH-1015
Lausanne
Switzerland 

\email{jerome.goudet@unil.ch}
}

\seealso{\code{\link{varcomp}}.}
\examples{
#load data set
data(gtrunchier)
attach(gtrunchier)
varcomp.glob(data.frame(Locality,Patch),gtrunchier[,-c(1,2)])
}
\keyword{univar}
