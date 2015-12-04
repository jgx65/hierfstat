\name{crocrussula}
\docType{data}
\alias{crocrussula}
\title{Genotypes and sex of 140 shrews Crocidura russula}
\description{A dataset containing microsatellite genotypes, population and sex of 140 Crocidura russula individuals}
\usage{
data(crocrussula)
}
\references{
 \href{http://rspb.royalsocietypublishing.org/content/264/1378/127}{Favre et al. (1997)} Female-biased dispersal in the monogamous mammal Crocidura russula: evidence from field data and microsatellite  patterns. Proceedings of the Royal Society, B (264):  127-132
 
\href{http://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.2002.01496.x/abstract}{Goudet J, Perrin N, Waser P (2002)} Tests for sex-biased dispersal using bi-parentally inherited genetic markers 11, 1103:1114
}

\examples{
data(crocrussula)
aic<-AIc(crocrussula$genot)
boxplot(aic~crocrussula$sex)
sexbias.test(crocrussula$genot,crocrussula$sex)
}

\keyword{datasets}

