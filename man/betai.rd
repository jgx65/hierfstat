\name{betai}
\alias{betai}
\title{Estimation of \eqn{\beta} per population}
\description{Estimates \eqn{\beta} (Fst) per population}
\usage{betai(dat)}
\arguments{
\item{dat}{A data frame containing the population of origin in the first column and the genotypes in the following ones}
}

\value{
\item{betai}{ \eqn{\beta_i} are from (7) of WH2002}
\item{betaio}{\eqn{\beta_{io}} are population estimates over loci 1-(sum of nums)/(sum of dens)}
\item{betaw}{\eqn{\beta_w} are average over populations using (9) of WH2002}
\item{betaii'}{\eqn{\beta_{ii'}} can be extracted with betas[1,,] etc...}
}
\details{These estimators assume alleles to be independant (i.e., random mating)}
\references{Weir and Hill (2002) ESTIMATING F-STATISTICS. Annual review of Genetics. 36:721}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
dat<-sim.genot(size=100,N=c(100,1000,10000),nbloc=50,nbal=10)
betai(dat)$betai 
}
\keyword{univar}
