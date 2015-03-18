\name{samp.within}
\alias{samp.within}
\title{Shuffles a sequence within groups defined by the input vector}
\description{Used to generate a permutation of a sequence \kbd{1:length(lev)}. observations are permutted within blocks, according to the vector \kbd{lev} passed to the function.}
\usage{
samp.within(lev)
}
\arguments{
\item{lev}{a vector containing the group to which belongs the observations to be permuted.}
}
\value{
a vector \kbd{1:length(lev)} (with blocks defined by \preformatted{lev}) randomly permuted. Usually, one passes the result to reorder observations in a data set in order to carry out permutation-based tests.
 }
\references{
Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186
}
\author{Jerome Goudet,
DEE, UNIL, 
CH-1015
Lausanne
Switzerland 

\email{jerome.goudet@unil.ch}
}
\seealso{\code{\link{samp.between}},\code{\link{g.stats.glob}}.}
\examples{
samp.within(rep(1:4,each=4))
#for an application see example in g.stats.glob
}
\keyword{manip}
