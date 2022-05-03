\name{getal.b}
\alias{getal.b}
\title{Converts diploid genotypic data into allelic data}
\description{Converts a data frame of genotypic diploid data with as many lines as
individuals (ni) and as many columns as loci (nl) into an array [ni,nl,2] of allelic data}
\usage{getal.b(data)}
\arguments{
\item{data}{a data frame with ni rows and nl columns. Each line encodes one individual, each column contains the genotype at one locus of the individual}
}

\value{
an array [ni,nl,2] of alleles. The two alleles are stored in the third dimension of the array
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
#multilocus diploid genotype of the first individual
gtrunchier[1,-c(1:2)]
#the diploid genotype splitted in its two constituent alleles
getal.b(gtrunchier[,-c(1:2)])[1,,]
}
\keyword{manip}
