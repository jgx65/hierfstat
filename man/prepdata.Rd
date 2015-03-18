\name{prepdata}
\alias{prepdata}
\title{Sort and renumber and recode levels for the hierarchical analysis if necessary}
\description{Called by varcomp to rearrange data. Should not need to be called directly}
\usage{
prepdata(data)
}
\arguments{
\item{data}{takes as input the different factor from the outermost (e.g. region) to the innermost (e.g. individual). This function must be applied before transformation of genotypic to allelic data}
}
\value{
a reoredered,renumbered and sorted data frame 
 }
\references{
Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and F-statistics. Molecular Ecology Notes. 5:184-186
}
\seealso{\code{\link{varcomp}}.}
\examples{
f1<-rep(c("B","C","A"),each=10)
f2<-rep(c("d","i","f","g","h","e"),each=5)
prepdata(data.frame(loc=f1,patch=f2,ind=1:30))
}
\keyword{manip}
