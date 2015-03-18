\name{nei.dist}
\alias{nei.dist}
\title{Estimates Nei's genetic distance}
\description{Estimates Nei's unbiased genetic distance}
\usage{nei.dist(data)}
\arguments{
\item{data}{a genetic data set, such as obtained from \kbd{read.fstat}}
}

\value{
\item{}{a vector of genetic distances among pairs 2-1, 3-1, 3-2,4-1,etc}
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{data(gtrunchier)
nei.dist(gtrunchier[,-1])}
\keyword{univar}
