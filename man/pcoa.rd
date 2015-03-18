\name{pcoa}
\alias{pcoa}
\title{Principal coordinate analysis}
\description{principal coordinates analysis as described in Legendre & Legendre Numerical Ecology}
\usage{pcoa(mat,plotit=TRUE,...)}
\arguments{
\item{mat}{a distance matrix}
\item{plotit}{Whether to produce a plot of the pcoa}
\item{...}{further arguments (graphical for instance) to pass to the function}
}

\value{
\item{valp}{the eigen values of the pcoa}
\item{vecp}{the eigen vectors of the pcoa (the coordinates of observations)}
\item{eucl}{The cumulative euclidian distances among observations, }
}
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{
data(gtrunchier)
colo<-c("black","red","blue","yellow","orange","green")
pcoa(vec2mat(cfe.dist(gtrunchier[,-1])),col=rep(colo,c(5,5,4,5,5,5)))
}

