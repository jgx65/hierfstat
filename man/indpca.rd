\name{indpca}
\alias{indpca}
\alias{print.indpca}
\alias{plot.indpca}
\title{PCA on a matrix of individuals genotypes frequencies}
\description{Carry out a PCA on the centered, unscaled matrix of individual's allele frequencies.}
\usage{
indpca(dat,ind.labels=NULL,scale=FALSE)

\method{print}{indpca}(x,...)
\method{plot}{indpca}(x,eigen=FALSE,ax1=1,ax2=2,...)
}
\arguments{
\item{dat}{A data frame with population of origin as first column, and genotypes in following columns.}
\item{ind.labels}{a vector of labels for the different individuals}
\item{scale}{whether to standardize each column to variance 1 or to leave it as is (default)}
\item{x}{an indpca object}
\item{eigen}{whether to plot in an additional windows screeplot of the inertias for the different axes}
\item{ax1}{which PCA coordinates to plot on the x axis}
\item{ax2}{which PCA coordinates to plot on the y axis}
\item{...}{further arguments to pass to print or plot}
}
\value{An object of class \code{indpca} with components
\item{call}{The function call}
\item{ipca}{an object of class pca and dudi (see dudi.pca) in package ade4}
\item{mati}{the original non centered matrice of individuals X alleles frequencies }
}

\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\examples{
##not run
data(gtrunchier)
x<-indpca(gtrunchier[,-2],ind.labels=gtrunchier[,2])
plot(x,col=gtrunchier[,1],cex=0.7)
}
\keyword{multivariate}
