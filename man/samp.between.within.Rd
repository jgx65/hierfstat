\name{samp.between.within}
\alias{samp.between.within}
\title{Shuffles a sequence}
\description{Used to generate a permutation of a sequence \kbd{1:length(inner.lev)}. blocks of observations defined by \kbd{inner.lev} are permutted within blocks defined by \kbd{outer.lev}}
\usage{
samp.between.within(inner.lev, outer.lev)
}
\arguments{
\item{inner.lev}{a vector containing the groups to be permuted.}
\item{outer.lev}{a vector containing teh blocks within which observations are to be kept.}
}
\value{
a vector \kbd{1:length(lev)} (with blocks defined by data) randomly permuted. Usually, one passes the result to reorder observations in a data set in order to carry out permutation-based tests
 }

\seealso{\code{\link{test.between.within}}.}
\keyword{manip}
