\name{basic.stats}
\alias{basic.stats}
\alias{print.basic.stats}
\title{Basic statistics}
\description{Estimates individual counts, allelic frequencies, observed heterozygosities and genetic diversities per locus and population.
Also Estimates mean observed heterozygosities, mean gene diversities within population Hs, Gene diversities overall Ht and corrected Htp, and Dst, Dstp.
Finally, estimates Fst and Fstp as well as Fis following Nei (1987) per locus and overall loci}
\usage{basic.stats(data,diploid=TRUE,digits=4)
\method{print}{basic.stats}(x,...)
}
\arguments{
\item{data}{a data frame where the first column contains the population to which the different individuals belong, and the following columns contain the genotype of the individuals -one locus per column- }
\item{diploid}{Whether individuals are diploids (default) or haploids}
\item{digits}{how many digits to print out in the output (default is 4)}
\item{x}{an object of class basic.stats}
\item{...}{further arguments to pass to print.bas.stats}
}

\value{
\item{n.ind.samp}{A table --with np (number of populations) columns and nl (number of loci) rows-- of genotype counts}
\item{pop.freq}{A list containing allele frequencies. Each element of the list is one locus.
For each locus, Populations are in columns and alleles in rows}
\item{Ho}{A table --with np (number of populations) columns and nl (number of loci) rows-- of observed heterozygosities}
\item{Hs}{A table --with np (number of populations) columns and nl (number of loci) rows-- of observed gene diversities}
\item{Fis}{A table --with np (number of populations) columns and nl (number of loci) rows--of observed Fis}
\item{perloc}{A table --with as many rows as loci-- containing basic statistics Ho, Hs, Ht, Dst, Ht', Dst', Fst, Fst' ,Fis, Dest}
\item{overall}{Basic statistics averaged over loci}
} 


\references{
Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press

Jost L (2008) GST and its relatives do not measure differentiation.
Molecular Ecology, 17, 4015-4026.

Nei M, Chesser R (1983) Estimation of fixation indexes and gene diversities.
Annals of Human Genetics, 47, 253-259.
}

\author{Jerome Goudet \email{jerome.goudet@unil.ch}}

\seealso{\code{\link{ind.count}},\code{\link{pop.freq}}.}

\note{
For the perloc and overall tables (see value section), the following statistics, defined in eq.7.38-- 7.43 pp.164--5 of Nei (1987)
are estimated:

The observed heterozygosity
\deqn{Ho= 1-\sum_k \sum_i  Pkii/np,}

where \eqn{Pkii} represents the proportion of homozygote \eqn{i} in sample \eqn{k} and
\eqn{np} the number of samples.

The within population gene diversity (sometimes misleadingly called expected
heterozygosity):

\deqn{Hs=\tilde{n}/(\tilde{n}-1)[1-\sum_i\bar{p_i^2}-Ho/2\tilde{n}],
}

where \eqn{\tilde{n}=np/\sum_k 1/n_k} and
\eqn{\bar{p_i^2}=\sum_k p_{ki}^2/np}

The overall gene diversity
\deqn{
Ht= 1-\sum_i\bar{p_i}^2+Hs/(\tilde{n} np)-Ho/(2\tilde{n}
np),}
where \eqn{\bar{p_i}=\sum_kp_{ki}/np}.

The amount of gene diversity among samples \eqn{Dst=Ht-Hs}

\eqn{Dst'=np/(np-1)Dst}

\eqn{Ht'=Hs+Dst'}

\eqn{Fst=Dst/Ht}.(This is not the same as Nei's \eqn{Gst},
Nei's \eqn{Gst} is an estimator of \eqn{Fst} based on allele frequencies only)

\eqn{Fst'=Dst'/Ht'}

\eqn{Fis=1-Ho/Hs}


Last, \eqn{Dest=np/(np-1) (Ht'-Hs)/(1-Hs)} a measure of population
differentiation as defined by Jost (2008) is also given

Here, the \eqn{p_{ki}} are unweighted by sample size. These statistics are
estimated for each locus and an overall loci estimates is also given, as the
unweighted average of the per locus estimates. In this way, monomorphic loci
are accounted for (with estimated value of 0) in the overall estimates.

Note that the equations used here all rely on genotypic rather than allelic
number and are corrected for heterozygosity.
}


\examples{
data(gtrunchier)
basic.stats(gtrunchier[,-1])
}
\keyword{univar}
