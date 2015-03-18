\name{eucl.dist}
\alias{eucl.dist}
\title{Estimates euclidian distances}
\description{Estimates euclidian distances among pairs of samples from a gnetic data set}
\usage{eucl.dist(data,allloc=FALSE,distance="eucl")}
\arguments{
\item{data}{a genetic data set, such as obtained from \kbd{read.fstat}}
\item{allloc}{whether to print out the distance for each locus}
\item{distance}{not used yet. Goal is to have only one distance function}
}

\value{
\item{}{either a matrix [np*(np-1)/2,nl] or  a vector [np*(np-1)/2] of genetic distances among pairs 2-1, 3-1, 3-2,4-1,etc}
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{data(gtrunchier)
nei.dist(gtrunchier[,-1])}
\keyword{univar}
