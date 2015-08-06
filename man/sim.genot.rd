\name{sim.genot}
\alias{sim.genot}
\title{Simulates genotypes in an island model at equilibrium}
\description{Simulates genotypes from several individuals in several populations at several loci in an island model at equilibrium. The islands may differ in size and inbreeding coeeficients.}
\usage{sim.genot(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0)}
\arguments{
\item{size}{The number of individuals to sample per population}
\item{nbal}{The maximum number of alleles present at a locus}
\item{nbloc}{The number of loci to simulate}
\item{nbpop}{The number of populations to simulate}
\item{N}{The population sizes for each island}
\item{mig}{the proportion of migration among islands}
\item{mut}{The loci mutation rate}
\item{f}{the inbreeding coefficient for each island}
}

\value{
a data frame with nbpop*size lines and nbloc+1 columns. Individuals are in rows and genotypes in columns, the first column being the population identifier}
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
dat<-sim.genot(nbpop=4,nbal=20,nbloc=10,mig=0.001,mut=0.0001,N=c(100,100,1000,1000),f=0)
betas(dat)$betaiovl
}
