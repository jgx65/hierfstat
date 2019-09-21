#' Converts a kinship matrix to a Genetic Relation Matrix (GRM)
#'  
#' @usage kinship2grm(x)
#' @param x a square matrix containing kinship coefficients
#' @return a GRM matrix
#' @details for off-diagonal elements, \eqn{GRM=2 \times x_{ij}}; for diagonal elements, \eqn{GRM=1+ x_{ii}}
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @examples
#' \dontrun{
#' dos<-matrix(sample(0:2,replace=TRUE,size=1000),nrow=10) #dosage matrix for 10 inds at 100 loci
#' ks<-beta.dosage(dos) # kinship matrix
#' kinship2grm(ks)
#' }  
#' @export 


kinship2grm<-function(x){
  fis<-diag(x)
  grm<-2*x
  diag(grm)<-1+fis
  return(grm)
}

#' Converts a kinship matrix to a distance matrix
#'  
#' @usage kinship2dist(x)
#' @param x A square matrix containg kinship coefficients
#' @value A distance matrix
#' @details \eqn{D_{ii}=0, D_{ij}=\frac{1-(x-min(x))}{(1-min(x))}}
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @export

kinship2dist<-function(x){
  diag(x)<-NA
  minx<-min(x,na.rm=TRUE)
  d.x<-1-(x-minx)/(1-minx)
  diag(d.x)<-0
  return(d.x)
}

#' Converts a Genetic Relationship Matrix (GRM) to a kinship matrix 
#'  
#' @usage grm2kinship(x)
#' @param x a square (GRM) matrix
#' @value a kinship matrix
#' @details \eqn{k[ii]=x[ii]-1; k[ij]=x[ij]/2}
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @export

grm2kinship<-function(x){
  fis<-diag(x)
  ks<-x/2
  diag(ks)<-fis-1
  return(ks)
}

#' Shifts a kinship matrix 
#' 
#' @usage kinshipShift(x,shift=NULL) 
#' @param x a square matrix
#' @param shift the amount by which the elements of x should be shifted. if \code{shift==NULL}, 
#' the average of the off-diagonal elements is substracted 
#' @value the shifted kinship matrix \eqn{\frac{x-shift}{1-shift}}
#' @details The kinship matrix produced by \code{beta.dosage} is relative to the average kinship
#' of the set of individuals analysed (\eqn{1/(n(n-1)/2) \sum_i \sum_{j>i} x_{ij}=0}). 
#' Another reference point might be useful, for instance to avoid negative kinship values, one might
#' want to shift the matrix by \eqn{min(x_{ij}), i \neq j}. 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @export

kinshipShift<-function(x,shift=NULL){
  if (is.null(shift)) shift<-mean(mat2vec(x))
  return((x-shift)/(1-shift))
}

#' Calculates population specific \eqn{F_{ST}} from dosage data  
#' 
#' @usage fst.dosage(dat,pop,matching=FALSE)
#' @param dat either a niXnl dosage matrix of 0,1 and 2s containing the 
#' number of derived alleles; or a niXni matrix of matching proportions
#' @param pop a vector of length ni containing the identifier of the population
#' to which individuals belong
#' @param matching a logical element, whether \code{dat} is dosage (niXnl, default) or 
#' matching (niXni)  matrix.
#' @value a vector of length np+1, the first np elements being population specific \eqn{F_{ST}}s, 
#' the last one being the overall \eqn{F_{ST}}
#' @note Does not assume Hardy-Weinberg equilibrium. If the numbers of samples per population 
#' are identical, fst.dosage gives the same results as \code{\link{wc}} for overall \eqn{F_{ST}}
#' @references 
#' 
#' \href{http://www.genetics.org/content/206/4/2085}{Weir, BS and Goudet J. 2017} A Unified Characterization 
#' of Population Structure and Relatedness. Genetics (2017) 206:2085 
#' 
#' @seealso The \code{\link{betas}} function estimates population specific \eqn{F_{ST}} from non-dosage data, 
#' assuming random mating
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @examples
#' \dontrun{ 
#' dos<-matrix(sample(0:2,replace=TRUE,size=1000),nrow=10) #dosage matrix for 10 inds at 100 loci
#' fst.dosage(dos,rep(1:2,each=5)) 
#' 
#' ks<-beta.dosage(dos,Mb=TRUE) # kinship matrix
#' fst.dosage(ks$betas*(1-ks$MB)+ks$MB,pop=rep(1:2,each=5),matching=TRUE)
#' 
#' dat<-sim.genot(size=100,nbal=2,nbloc=200)
#' c(betas(dat)$betaiovl,betas(dat)$betaW)
#' dos<-as.matrix((dat[,-1]%/%10-1)+(dat[,-1]%%10-1))
#' fst.dosage(dos,pop=rep(1:3,each=100))
#' wc(dat)$FST # same as fst.dosage overall if sampling is balanced 
#' }
#' @export

fst.dosage<-function(dat,pop,matching=FALSE){
  if(!matching){
    nl <- dim(dat)[2]
    Mij <- 1/2 + tcrossprod(dat - 1)/2/nl
    Mii <- (diag(Mij)) * 2 - 1
    diag(Mij) <- NA
  }
  else {
    Mij<-dat
    Mii<-diag(Mij)
    diag(Mij)<-NA
  }
  pop<-factor(pop)
  x<-levels(pop)
  wil<-lapply(x,function(z) which(pop==z))
  Fsts<-unlist(lapply(wil,function(x) mean(Mij[x,x],na.rm=T)))
  
  for (i in 1:length(wil)){
    z<-wil[[i]]
    Mij[z,z]<-NA
  }
  Mb<-mean(Mij,na.rm=T)
  Fsts<-(Fsts-Mb)/(1-Mb)
  res<-c(Fsts,mean(Fsts,na.rm=TRUE))
  names(res)<-c(levels(pop),"All")
  return(res)
}


