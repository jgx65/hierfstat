#' @title Estimates pairwise kinships and individual inbreeding coefficients from dosage data
#'
#' @description Estimates pairwise kinships (coancestries) and individual inbreeding coefficient 
#' using Weir and Goudet (2017) beta estimator. 
#'  
#' @usage beta.dosage(dat,inb=TRUE,correction=FALSE,Mb=FALSE)
#'
#' @param dat A matrix of 0, 1 and 2s with loci (SNPs) in columns and individuals in rows. missing values are allowed
#' @param inb whether individual inbreeding coefficient should be estimated
#' @param correction whether kinships and inbreeding coefficients should be corrected following Goudet et al. (2018)
#' @param Mb whether to output the mean matching 
#' 
#' @return either a matrix of pairwise kinships and inbreeding coefficients if requested (if \code{Mb}=FALSE) or
#' a list with elements \code{inb} (whether inbeeding coefficients rather than kinships should be returned on the main diagonal),
#' \code{correction} (whether to return kinship or inbreeding corrected for missing values),
#' \code{MB} (the average matching) and \code{betas} the kinships or inbreeding coefficients.
#'
#' @details This function is written for dosage data, i.e., how many doses of an allele (0, 1 or 2) an individual carries.
#' It should be use for bi-allelic markers only (e.g. SNPs), although you might "force" a k multiallelic locus to k biallelic
#' loci. 
#' 
#' Matching probabilities can be obtained by the following equation: \eqn{M=\beta*(1-Mb)+Mb} 
#' 
#' When there are missing data, the missing values are replaced by the frequency of the alternate allele 
#' The correction option unbiases the estimates, and is described in the supplementary materials of Goudet etal. (2018). 
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references 
#' 
#' \href{http://www.genetics.org/content/206/4/2085}{Weir, BS and Goudet J. 2017} A Unified Characterization 
#' of Population Structure and Relatedness. Genetics (2017) 206:2085 
#' 
#' \href{https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14833}{Goudet, J., Kay, T. and Weir BS. 2018} How to estimate kinship. 
#' Molecular Ecology 27:4121.
#'
#' @examples 
#' \dontrun{
#'  dos<-matrix(sample(0:2,size=10000,replace=TRUE),ncol=100)
#'  beta.dosage(dos,inb=TRUE)
#' }
#' @export

beta.dosage<-function(dat,inb=TRUE,correction=FALSE,Mb=FALSE){
  #dat is a data frame with individuals in rows and allelic dosage for each locus in colums  
  #uses matching proba -same equation as for population i.e. Mij=[xiXj+(2-xi)(2-xj)]/4
  #missing values are replaced by mean allelic dosage of the locus
  dat<-as.matrix(dat)
  lims<-range(dat,na.rm=TRUE)
  if ((lims[2]>2) | (lims[1]<0)) stop("input dosage matrix should contains only 0, 1 and 2s")
  
  prop.miss<-sum(is.na(dat))/prod(dim(dat))  
  if (prop.miss>0.0){
    i.miss<-rowMeans(is.na(dat)) 
    cm<-colMeans(dat,na.rm=TRUE)
    dum<-sweep(dat,2,cm,FUN="-")
    dum[is.na(dum)]<-0.0
    dat<-sweep(dum,2,cm,FUN="+")
  }  else i.miss<-numeric(nrow(dat))

  
  nl<-dim(dat)[2]
  Mij<-1/2+tcrossprod(dat-1)/2/nl 
  Mii<-diag(Mij)
  diag(Mij)<-NA
  MB<-mean(Mij,na.rm=T)
  diag(Mij)<-Mii
  coa<-(Mij-MB)/(1-MB)
  if (inb) {
    ind.inb<-((Mii*2-1)-MB)/(1-MB)
    diag(coa)<-ind.inb
  }
  if (correction) {
    correc<-outer((1-i.miss),(1-i.miss))
    ic<-diag(coa)
    coa<-coa/correc 
    if (inb) {
    ic<-ic + i.miss
    diag(coa)<-ic*(1.0-i.miss)
    }
    }  
  if (!Mb) return(coa) else return(list(inb=inb,correction=correction,MB=MB,betas=coa))

}
