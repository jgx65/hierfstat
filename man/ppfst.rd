\name{pp.fst}
\alias{pp.fst}
\title{fst per pair}
\description{fst per pair following Weir and Cockerham (1984)}
\usage{pp.fst(dat=dat,diploid=TRUE,...)}
\arguments{
\item{dat}{a genetic data frame}
\item{diploid}{whether data from diploid organism}
\item{...}{further arguments to pass to the function}
}

\value{
\item{call}{function call}
\item{fst.pp}{pairwise Fsts}
\item{vc.per.loc}{for each pair of population, the variance components per locus}
}
\references{
Weir, B.S. and Cockerham C.C. (1984) Estimating F-Statistics for the Analysis of Population Structure 38:1358-1370 \url{http://www.jstor.org/stable/2408641}

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.
}

\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
%\examples{}
\keyword{univar}
