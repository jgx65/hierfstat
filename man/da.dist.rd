\name{da.dist}
\alias{da.dist}
\title{Estimates Nei's DA distance}
\description{Estimates Nei's DA distance}
\usage{da.dist(data,allloc=FALSE,distance="da")}
\arguments{
\item{data}{a gnetic data set, such as obtained from \kbd{read.fstat}}
\item{allloc}{whether to print out the distance for each locus}
\item{distance}{not used yet. Goal is to have only one distance function}
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
