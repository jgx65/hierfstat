\name{write.struct}
\alias{write.struct}
\title{Write structure file}
\description{Write a genotype data set to a file in the structure format}
\usage{write.struct(dat,ilab=NULL,pop=NULL,MARKERNAMES=FALSE,MISSING=-9,fname="dat.str")}
\arguments{
\item{dat}{a genotype dataframe}
\item{ilab}{an (optional) column with individual labels}
\item{pop}{an (optional) column with population identifiers}
\item{MARKERNAMES}{whether to add a row  with marker names. If TRUE, takes the loci names from dat}
\item{MISSING}{The code for missing alleles}
\item{fname}{a string containing the file name (default to "dat.str")}}
\value{a text file in the structure format}
\references{Pritchard JK etal. 2000. Inference of population structure using multilocus genotype data. Genetics 155:945-959}
\author{Jerome Goudet \email{jerome.goudet@unil.ch}}
\keyword{miscellaneous}