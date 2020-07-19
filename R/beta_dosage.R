#' Estimates pairwise kinships and individual inbreeding coefficients from dosage data
#'
#' Estimates pairwise kinships (coancestries) and individual inbreeding coefficient 
#' using Weir and Goudet (2017) beta estimator. 
#' 
#' This function is written for dosage data, i.e., how many doses of an allele (0, 1 or 2) an individual carries.
#' It should be use for bi-allelic markers only (e.g. SNPs), although you might "force" a k multiallelic locus to k biallelic
#' loci. 
#' 
#' Matching proportion can be obtained by the following equation: \eqn{M=\beta*(1-Mb)+Mb} 
#' 
#' By default (inb=TRUE) the inbreeding coefficient is returned on the main diagonal.  With inb=FALSE, self coancestries are reported. 
#' 
#' @usage beta.dosage(dos,inb=TRUE,Mb=FALSE)
#'
#' @param dos A matrix of 0, 1 and 2s with loci (SNPs) in columns and individuals in rows. missing values are allowed
#' @param inb whether individual inbreeding coefficient should be estimated (rather than self-coancestries)
#' @param Mb whether to output the mean matching 
#' 
#' @return if \code{Mb}=FALSE, a matrix of pairwise kinships and inbreeding coefficients (if \code{inb}=TRUE) or self-coancestries 
#' (\code{inb}=FALSE);  
#' if \code{Mb}=FALSE, a list with elements \code{inb} (whether inbreeding coefficients rather than kinships should 
#' be returned on the main diagonal),
#' \code{MB} (the average off-diagonal matching) and \code{betas} the kinships or inbreeding coefficients.
#'
#' @details 
#' Twice the betas with self-coancestries on the diagonal gives the Genomic Relationship Matrix (GRM)     
#' 
#' Following a suggestion from Olivier Hardy, missing data are removed from the estimation procedure, rather than imputed (this is taken care off automatically) 
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


beta.dosage <- function (dos, inb = TRUE, Mb = FALSE) {
  #dos is a data frame with individuals in rows and allelic dosage for each locus in colums  
  #uses matching proba -same equation as for population i.e. Mij=[xiXj+(2-xi)(2-xj)]/4
  if(!is.matrix(dos)){
  if(class(dos)[[1]]=="bed.matrix") dos<-gaston::as.matrix(dos) 
  else dos <- as.matrix(dos)
  }
  lims <- range(dos, na.rm = TRUE)
  if ((lims[2] > 2) | (lims[1] < 0)) 
    stop("input dosage matrix should contains only 0, 1 and 2s")
  
  na <- matrix(rep(1,prod(dim(dos))),ncol=ncol(dos))
  ina<-which(is.na(dos))
  na[ina]<-0
  dos[ina]<-1
  
  Mij <- 1/2 * (1+1/tcrossprod(na) * tcrossprod(dos-1)) 
  
  Mii <- diag(Mij)
  diag(Mij) <- NA
  MB <- mean(Mij, na.rm = T)
  diag(Mij) <- Mii
  coa <- (Mij - MB)/(1 - MB)
  if (inb) {
    ind.inb <- ((Mii * 2 - 1) - MB)/(1 - MB)
    diag(coa) <- ind.inb
  }
  if (!Mb) 
    return(coa)
  else return(list(inb = inb, MB = MB, 
                   betas = coa))
}


