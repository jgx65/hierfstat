\name{read.fstat}
\alias{read.fstat}
\title{Reads data from a FSTAT file}
\description{Imports a \emph{FSTAT} data file into R.
The data frame created is made of \kbd{nl+1} columns, \kbd{nl} being the number of loci.
The first column corresponds to the Population identifier,
the following columns contains the genotypes of the individuals.}
\usage{
read.fstat(fname, na.s = c("0","00","000","0000","00000","000000","NA"))
}
\arguments{
\item{fname}{a file in the FSTAT format (\url{http://www.unil.ch/popgen/softwares/fstat.htm}):
The file must have the following format:

The first line contains 4 numbers: the number of samples, \kbd{np} , the number of loci, \kbd{nl},
the highest number used to label an allele, \kbd{nu}, and a 1 if the code for alleles is a one digit
 number (1-9), a 2 if code for alleles is a 2 digit number (01-99) or a 3 if code for alleles is a 3
 digit number (001-999).  These 4 numbers need to be separated by any number of spaces.}

The first line is immediately followed by \kbd{nl} lines, each containing the name of a locus,
in the order they will appear in the rest of the file.

	On line \kbd{nl+2}, a series of numbers as follow:

	\preformatted{1     0102   0103   0101  0203          0      0303}

The first number identifies the sample to which the individual belongs,
the second is the genotype of the individual at the first locus,
coded with a 2 digits number for each allele,
the third is the genotype at the second locus, until locus nl is entered
(in the example above, \kbd{nl=6}).  Missing genotypes are encoded with 0, 00, 0000, 000000 or NA.
Note that 0001 or 0100 are not a valid format, as both alleles at a locus have
to be known, otherwise, the genotype is considered as missing.
No empty lines are needed between samples.
\item{na.s}{The strings that correspond to the missing value.
\emph{You should note have to change this}}
}
\value{
a data frame containing the desired data, in a format adequate to pass to \kbd{varcomp}
 }
\references{
Goudet J. (1995). FSTAT (Version 1.2): A computer program to calculate F- statistics. Journal of Heredity 86:485-486

Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186
}
\examples{
read.fstat(paste(path.package("hierfstat"),"/extdata/diploid.dat",sep="",collapse=""))
}

\keyword{manip}
