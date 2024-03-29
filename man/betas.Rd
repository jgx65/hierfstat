% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betas.R
\name{betas}
\alias{betas}
\alias{print.betas}
\title{Estimates \eqn{\beta}s per population and a bootstrap confidence interval}
\usage{
betas(dat,nboot=0,lim=c(0.025,0.975),diploid=TRUE,betaijT=FALSE)

\method{print}{betas}(x, digits = 4, ...)
}
\arguments{
\item{dat}{data frame with genetic data and pop identifier}

\item{nboot}{number of bootstrap samples.}

\item{lim}{width of the bootstrap confidence interval}

\item{diploid}{whether the data comes from a diploid organism}

\item{betaijT}{whether to estimate individual coancestries}

\item{x}{a betas object}

\item{digits}{number of digits to print}

\item{...}{further arguments to pass to print}
}
\value{
Hi Within population gene diversities (complement to 1 of matching probabilities)

Hb Between populations gene diversities

betaiovl Average \eqn{\hat{\beta_{WT}^i}} over loci (Population specific FSTs), Table 3 of 
\href{https://academic.oup.com/genetics/article/206/4/2085/6072590}{Weir and Goudet, 2017 (Genetics)}

betaW Average of the betaiovl \eqn{\hat{\beta_{WT}}} over loci (overall population FST)

ci The bootstrap confidence interval of population specific FSTs
(only if more than 100 bootstraps requested AND if more than 10 loci are present)

if betaijT=TRUE, return the matrix of pairwise kinships only.
}
\description{
Estimate populations (Population specific FST) or individual coancestries 
and a bootstrap confidence interval, assuming random mating
}
\details{
If betaijT=TRUE, and the first column contains a unique identifier for 
each individual, the function returns the matrix of individual coancestries/kinships.  
Individual inbreeding coefficients can be obtained by multiplying by 2 the diagonal 
and substracting 1.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(betas)}: print function for betas class

}}
\examples{
\dontrun{
#3 different population sizes lead to 3 different betais
dat<-sim.genot(size=40,N=c(50,200,1000),nbloc=50,nbal=10)
betas(dat,nboot=100)
 
#individual coancestries from the smallest population are large
ind.coan<-betas(cbind(1:120,dat[,-1]),betaij=T)
diag(ind.coan$betaij)<-NA
graphics::image(1:120,1:120,ind.coan$betaij,xlab="Inds",ylab="Inds")
}

}
\references{
\href{https://academic.oup.com/genetics/article/206/4/2085/6072590}{Weir and Goudet, 2017 (Genetics)} 
A unified characterization of population structure and relatedness.
}
\seealso{
\code{\link{fs.dosage}}, \code{\link{beta.dosage}} for Fst estimates (not assuming Random Mating) 
and kinship estimates from dosage data, respectively
}
\author{
Jerome Goudet \email{jerome.goudet@unil.ch}
}
