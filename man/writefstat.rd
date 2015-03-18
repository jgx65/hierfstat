\name{write.fstat}
\alias{write.fstat}
\title{Write an Fstat data file}
\description{Write a data frame to a text file in the fstat data format, see \code{read.fstat}}
\usage{write.fstat(dat,fname="genotypes.dat")}
\arguments{
\item{dat}{A data frame with first column containing the population identifier and remaining columns containing genotypes}
\item{fname}{The name of teh text file to which the data frame should be written}
}
\value{None}
\references{
Goudet J. (1995). FSTAT (Version 1.2): A computer program to calculate F- statistics. Journal of Heredity 86:485-486
}
\author{Jerome Goudet}
\examples{
\dontrun{data(gtrunchier)
write.fstat(gtrunchier[,-1],"galba.dat")
}}
\keyword{IO}