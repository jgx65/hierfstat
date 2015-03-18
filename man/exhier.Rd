\name{exhier}
\docType{data}
\alias{exhier}
\title{Example data set with 4 levels, one diploid and one haploid locus}
\description{Example data set with 4 levels, one diploid and one haploid locus}
\usage{
data(exhier)
}
\value{
\item{lev1}{outermost level}
\item{lev2}{level 2}
\item{lev3}{Level 3}
\item{lev4}{Level 4}
\item{diplo}{Diploid locus}
\item{haplo}{Haploid locus}

}

\examples{
data(exhier)
varcomp(exhier[,1:5])
varcomp(exhier[,c(1:4,6)],diploid=FALSE)

}
\keyword{datasets}
