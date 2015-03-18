\name{getal}
\alias{getal}
\title{Converts diploid genotypic data into allelic data}
\description{Converts diploid genotypic data into allelic data}
\usage{getal(data)}
\arguments{
\item{data}{a data frame where the first column contains the population to which the different individuals belong, and the following columns contain the genotype of the individuals -one locus per column- }
}

\value{
\item{data.al}{a new data frame, with twice as many row as the input data frame and one extra column. each row of the first half of the data frame contains the first allele for each locus, and each row of the second half of the data frame contains the second allel at the locus. The extra column in second position corresponds to the identifier of the individual to which the allele belongs}
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
getal(data.frame(gtrunchier[,-2]))
}
\keyword{manip}
