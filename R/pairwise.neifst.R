###################################################
#'
#' @title Estimate pairwise FSTs according to Nei (1987)
#'
#' @description Estimate pairwise FSTs according to Nei (1987)
#' 
#' @usage pairwise.neifst(dat,diploid=TRUE)
#' 
#' @param dat A data frame containing population of origin as the first column and multi-locus genotypes in following columns
#' @param diploid whether the data is from a diploid (default) or haploid organism
#' 
#' @return A matrix of pairwise FSTs
#' 
#' @details FST are calculated using Nei (87) equations for FST', as described in the note section of  \link{basic.stats}
#' 
#' @references Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University Press
#' 
#'           
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @seealso \link{pp.fst} \link{genet.dist}
#' 
#' 
#' @examples
#' 
#' data(gtrunchier)
#' pairwise.neifst(gtrunchier[,-2],diploid=TRUE)
#' 
#' @export
#' 
#########################################################################

pairwise.neifst <- function(dat,diploid=TRUE){
  dat<-dat[order(dat[,1]),]
  pops<-unique(dat[,1])
  npop<-length(pops)
  fstmat <- matrix(nrow=npop,ncol=npop,dimnames=list(pops,pops))
  for(a in pops[-1]){
    for(b in pops[1]:(a-1)){ #pb if levels of a are factors
      subdat <- dat[dat[,1] %in% c(a,b),]
      fstmat[a,b]<-fstmat[b,a]<- basic.stats(subdat,diploid=diploid)$overall[8]
    }
    
  }
fstmat
}
