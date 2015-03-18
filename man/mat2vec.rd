\name{mat2vec}
\alias{mat2vec}
\title{Rewrite a matrix as a vecor}
\description{transform lower triangular matrix in vector 1.2,1.3,2.3,1.4,2.4,3.4}
\usage{mat2vec(mat)}
\arguments{
\item{mat}{a (square) matrix}
}

\value{
\item{x}{a vector}
} 
%\references{}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

%\seealso{\code{\link{}}.}
\examples{mat2vec(matrix(1:16,nrow=4))}
