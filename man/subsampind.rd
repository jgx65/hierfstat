\name{subsampind}
\alias{subsampind}
\title{Subsample a FSTAT data frame}
\description{Subsample a given number of individuals from a FSTAT data frame}
\usage{subsampind(dat,sampsize = 10)}
\arguments{
\item{dat}{A data frame with population of origin as first column, and genotypes in following columns.}
\item{sampsize}{the number of individuals to sample in each population.}
}
\value{A data frame with population of origin as first column, and genotypes in following columns.  Each population is made of at most sampsize individuals
}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\examples{
data(gtrunchier)
subsampind(gtrunchier[,-1],6)  # check the warning
}
\keyword{utilities}