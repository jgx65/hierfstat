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
#' @return A distance matrix
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
#' @return a kinship matrix
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
#' @return the shifted kinship matrix \eqn{\frac{x-shift}{1-shift}}
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


