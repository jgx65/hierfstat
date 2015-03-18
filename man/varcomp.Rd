\name{varcomp}
\alias{varcomp}
\title{Estimates variance components for each allele of a locus}
\description{Estimates variance components for each allele for a (fully) hierarchical random design defined by all but the last column of the data frame  \kbd{data}, the last column containing the genetic data to analyse. Columns for the hierarchical design should be given from the outermost to the innermost before the individual (e.g. continent, region, population, patch,...)}
\usage{
varcomp(data,diploid=TRUE)
}
\arguments{
\item{data}{a data frame that contains the different factors from the outermost (e.g. region)
to the innermost before the individual. the last column of the data frame 'data' contains the locus to analyse, which can be multiallelic. Missing data are allowed.}
\item{diploid}{a boolean stating whether the data come from diploid (TRUE=default) or haploid (FALSE) organisms}}
\value{
\item{df}{the degrees of freedom for each level}
\item{k}{the k matrix, the coefficients associated with the variance components}
\item{res}{the variance components for each allele}
\item{overall}{the variance components summed over alleles}
\item{F}{a matrix of hierarchical F-statistics type-coefficients
 with the first line corresponding to \eqn{F_{(n-1)/n}},\eqn{F_{(n-2)/n}}...\eqn{F_{i/n}} and
 the diagonal corresponding to \eqn{F_{(n-1)/n}}, \eqn{F_{(n-2)/(n-1)}},\eqn{F_{i/2}}}
}
\details{The format for genotypes is simply the code for the 2 alleles put one behind the other, without space in between. For instance if allele 1 at the locus has code \kbd{23} and allele 2 \kbd{39}, the genotype format is \kbd{2339}. 
}
\references{
Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52(4):950-956
}
\author{Jerome Goudet,
DEE, UNIL, 
CH-1015
Lausanne
Switzerland 

\email{jerome.goudet@unil.ch}

\url{http://www.unil.ch/popgen/people/jerome.htm}
}
\seealso{\code{\link{varcomp.glob}}.}
\examples{
#load data set
data(gtrunchier)
attach(gtrunchier)
#
varcomp(data.frame(Locality,Patch,L21.V))
}
\keyword{univar}
