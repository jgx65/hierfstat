\name{boot.ppfst}
\alias{boot.ppfst}
\alias{print.boot.ppfst}
\title{Performs bootstrapping over loci of pairwise Fst}
\description{Performs bootstrapping over loci of pairwise Fst using Weir and Cockerham (1984) estimator of Fst}
\usage{boot.ppfst(dat=dat,nboot=100,quant=c(0.025,0.975),diploid=TRUE,...)}
\arguments{
\item{dat}{a genetic data frame with population identifiers in the first
column. Numeric, character and factor identifiers are supported. At least two
observed populations are required; identifiers must not be missing. Unused
factor levels are ignored.}
\item{nboot}{number of bootstraps}
\item{quant}{the quantiles for bootstrapped ci}
\item{diploid}{whether data are from diploid organisms}
\item{...}{further arguments to pass to the function}
}
\value{
\item{call}{call to the function}
\item{ll}{lower limit ci}
\item{ul}{upper limit ci}
\item{vc.per.loc}{for each pair of population, the variance components per locus}
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
x<-boot.ppfst(gtrunchier[,-2])
x$ll
x$ul
}
