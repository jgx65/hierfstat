\name{g.stats.glob}
\alias{g.stats.glob}
\title{Likelihood ratio G-statistic over loci}
\description{Calculates the likelihood ratio G-statistic on a contingency table of alleles at one locus X sampling unit, and sums this statistic over the loci provided. The sampling unit could be any hierarchical level (patch, locality, region,...). By default, diploid data are assumed}
\usage{
g.stats.glob(data,diploid=TRUE)
}
\arguments{
\item{data}{a data frame made of \kbd{nl+1} column, \kbd{nl} being the number of loci. The first column contains the sampling unit, the others the multi-locus genotype. Only complete multi-locus genotypes are kept for calculation}
\item{diploid}{Whether the data are from diploid (default) organisms}
}
\value{
\item{g.stats.l}{Per locus likelihood ratio statistic}
\item{g.stats}{Overall loci likelihood ratio statistic}
 }
\references{
Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186

Goudet J., Raymond, M., DeMeeus, T. and Rousset F. (1996) Testing differentiation in diploid populations. Genetics. 144: 1933-1940

Petit E., Balloux F. and Goudet J.(2001) Sex-biased dispersal in a migratory bat: A characterization using sex-specific demographic parameters. Evolution 55: 635-640. 
}
\author{Jerome Goudet,
DEE, UNIL, 
CH-1015
Lausanne
Switzerland 

\email{jerome.goudet@unil.ch}
}
\seealso{\code{\link{g.stats}}, \code{\link{samp.within}},\code{\link{samp.between}}.}
\examples{
\dontrun{
data(gtrunchier)
attach(gtrunchier)
nperm<-99
nobs<-length(Patch)
gglobs.o<-vector(length=(nperm+1))
gglobs.p<-vector(length=(nperm+1))
gglobs.l<-vector(length=(nperm+1))

gglobs.o[nperm+1]<-g.stats.glob(data.frame(Patch,gtrunchier[,-c(1,2)]))$g.stats
gglobs.p[nperm+1]<-g.stats.glob(data.frame(Patch,gtrunchier[,-c(1,2)]))$g.stats
gglobs.l[nperm+1]<-g.stats.glob(data.frame(Locality,gtrunchier[,-c(1,2)]))$g.stats

for (i in 1:nperm) #careful, might take a while
{
  gglobs.o[i]<-g.stats.glob(data.frame(Patch,gtrunchier[sample(Patch),-c(1,2)]))$g.stats
  gglobs.p[i]<-g.stats.glob(data.frame(Patch,gtrunchier[samp.within(Locality),-c(1,2)]))$g.stats
  gglobs.l[i]<-g.stats.glob(data.frame(Locality,gtrunchier[samp.between(Patch),-c(1,2)]))$g.stats
}
#p-value of first test (among patches)
p.globs.o<-sum(gglobs.o>=gglobs.o[nperm+1])/(nperm+1) 

#p-value of second test (among patches within localities)
p.globs.p<-sum(gglobs.p>=gglobs.p[nperm+1])/(nperm+1) 

#p-value of third test (among localities)
p.globs.l<-sum(gglobs.l>=gglobs.l[nperm+1])/(nperm+1) 


#Are alleles associated at random among patches
p.globs.o 

#Are alleles associated at random among patches within localities?
#Tests differentiation among patches within localities
p.globs.p 

#Are alleles associated at random among localities, keeping patches as one unit?
#Tests differentiation among localities
p.globs.l 
}
}
\keyword{univar}
