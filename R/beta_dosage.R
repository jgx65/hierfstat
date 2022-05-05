#' Estimates matching between pairs of individuals
#' 
#' Estimates matching between pairs of individuals (for each locus, gives 1 if the two individuals are homozygous
#' for the same allele, 0 if they are homozygous for a different allele, and 1/2 if at least one individual 
#' is heterozygous. Matching is the average of these 0, 1/2 and 1s)
#'
#' This function is written for dosage data, i.e., how many doses of an allele (0, 1 or 2) an individual carries.
#' It should be use for bi-allelic markers only (e.g. SNPs), although you might "force" a k multiallelic locus to k 
#' biallelic loci (see \code{\link{fstat2dos}}). 
#'
#' @usage matching(dos)
#'
#' @param dos A matrix of 0, 1 and 2s with loci (SNPs) in columns and individuals in rows. 
#' missing values are allowed
#' 
#' @return a matrix of pairwise matching
#' 
#' @export

matching<-function(dos){
  if(!is.matrix(dos)){
    if(class(dos)[[1]]=="bed.matrix") dos<-gaston::as.matrix(dos) 
    else dos <- as.matrix(dos)
  }
  lims <- range(dos, na.rm = TRUE)
  if ((lims[2] > 2) | (lims[1] < 0)) 
    stop("input dosage matrix should contains only 0, 1 and 2s")
  
  if(sum(is.na(dos))>0){
    na <- matrix(rep(1,prod(dim(dos))),ncol=ncol(dos))
    ina<-which(is.na(dos))
    na[ina]<-0
    dos[ina]<-1
    Mij <- 1/2 * (1+1/tcrossprod(na) * tcrossprod(dos-1)) 
  }
  else {
    nl<-dim(dos)[2]
    Mij<-1/2 * (1+tcrossprod(dos-1)/nl)
  }
Mij  
}
  


#' Estimates pairwise kinships and individual inbreeding coefficients from dosage data
#'
#' Estimates pairwise kinships (coancestries) and individual inbreeding coefficient 
#' using Weir and Goudet (2017) beta estimator. 
#' 
#' This function is written for dosage data, i.e., how many doses of an allele (0, 1 or 2) an individual carries.
#' It should be use for bi-allelic markers only (e.g. SNPs), although you might "force" a k multiallelic locus to k biallelic
#' loci (see \code{\link{fstat2dos}}). 
#' 
#' Matching proportions can be obtained by the following equation: \eqn{M=\beta*(1-Mb)+Mb} 
#' 
#' By default (\code{inb=TRUE}) the inbreeding coefficient is returned on the main diagonal.  
#' With \code{inb=FALSE}, self coancestries are reported. 
#' 
#' @usage beta.dosage(dos,inb=TRUE,Mb=FALSE,MATCHING=FALSE)
#'
#' @param dos A matrix of 0, 1 and 2s with loci (SNPs) in columns and individuals in rows. 
#' Missing values are allowed
#' @param inb whether individual inbreeding coefficient should be estimated (rather than self-coancestries)
#' @param Mb whether to output the mean matching 
#' @param MATCHING if \code{MATCHING=FALSE}, \code{dos} is a (ni x nl) dosage matrix; 
#' if \code{MATCHING=TRUE}, dos is a (ni x ni) matrix 
#' of matching proportions, as obtained from a call to the \code{\link{matching}} function 
#' 
#' @return if \code{Mb}=FALSE, a matrix of pairwise kinships and inbreeding coefficients (if \code{inb}=TRUE) or self-coancestries 
#' (\code{inb}=FALSE);  
#' if \code{Mb}=TRUE, a list with elements \code{inb} (whether inbreeding coefficients rather than kinships should 
#' be returned on the main diagonal),
#' \code{MB} (the average off-diagonal matching) and \code{betas} the kinships or inbreeding coefficients.
#'
#' @details 
#' Twice the betas with self-coancestries on the diagonal gives the Genomic Relationship Matrix (GRM)     
#' 
#' Following a suggestion from Olivier Hardy, missing data are removed from the estimation procedure, 
#' rather than imputed (this is taken care off automatically) 
#'
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @references 
#' 
#' \href{https://academic.oup.com/genetics/article/206/4/2085/6072590}{Weir, BS and Goudet J. 2017} A Unified Characterization 
#' of Population Structure and Relatedness. Genetics (2017) 206:2085 
#' 
#' \href{https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14833}{Goudet, J., Kay, T. and Weir BS. 2018} How to estimate kinship. 
#' Molecular Ecology 27:4121.
#'
#' @examples 
#' \dontrun{
#'  dos<-matrix(sample(0:2,size=10000,replace=TRUE),ncol=100)
#'  beta.dosage(dos,inb=TRUE)
#'  
#'  #matrix of kinship/inbreeding coeff
#'  data(gtrunchier)
#'  beta.dosage(fstat2dos(gtrunchier[,-c(1:2)]))
#'  
#'  #individual inbreeding coefficients
#'  dat<-sim.genot(size=100,nbloc=100,nbal=20,mig=0.01,f=c(0,0.3,0.7))
#'  hist(diag(beta.dosage(fstat2dos(dat[,-1]))),breaks=-10:100/100,main="",xlab="",ylab="")
#'  abline(v=c(0.0,0.3,0.7),col="red")
#'  #only 20 loci
#'  hist(diag(beta.dosage(fstat2dos(dat[,2:21]))),breaks=-5:20/20,main="",xlab="",ylab="")
#'  abline(v=c(0.0,0.3,0.7),col="red")

#'  
#' }
#' @export


beta.dosage <- function (dos, inb = TRUE, Mb = FALSE, MATCHING=FALSE) {
  #dos is a data frame with individuals in rows and allelic dosage for each locus in colums  
  #uses matching proba -same equation as for population i.e. Mij=[xiXj+(2-xi)(2-xj)]/4

  if(!MATCHING)  {Mij<-matching(dos)} else Mij<-dos
  

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


