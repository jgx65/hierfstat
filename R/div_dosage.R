#' @title Estimates nucleotide diversity (\eqn{\pi}) from dosage data
#' @description Estimates nucleotide diversity \eqn{\pi= \sum_l 2 p_ l(1-p_l) 2n/(2n-1)} from a
#' dosage matrix
#' @usage pi.dosage(dos,L=NULL)
#' @param dos a ni X nl dosage matrix containing the number of derived/alternate alleles each individual carries
#' at each SNP 
#' @param L the length of the sequence
#' @return if \code{L=NULL} (default), returns the sum over SNPs of nucleotide diversity; 
#' otherwise return the average nucleotide diversity per nucleotide given the length \code{L} of the sequence  
#' 
################################################### 
#' @export

pi.dosage<-function(dos,L=NULL){
  if(!is.matrix(dos)){
    if(class(dos)[[1]]=="bed.matrix") dos<-gaston::as.matrix(dos) 
    else dos <- as.matrix(dos)
  }
  lims <- range(dos, na.rm = TRUE)
  if ((lims[2] > 2) | (lims[1] < 0)) 
    stop("input dosage matrix should contains only 0, 1 and 2s")

nis<-nrow(dos)  
p<-colMeans(dos,na.rm=TRUE)/2
nas<-colSums(is.na(dos))
pis<-2*p*(1-p)*(2*(nis-nas))/(2*(nis-nas)-1)
if(is.null(L)) sum(pis,na.rm=TRUE) else sum(pis,na.rm=TRUE)/L

}


#####################################################
#' @title Estimates \eqn{\theta_{Watterson}} from dosage data
#' @description Estimates \eqn{\theta_{Watterson}=S/a}, where \eqn{S} is the number of segregating sites
#' in a set of sequences and \eqn{a=1/\sum_i^{n-1} i}.  
#' @usage theta.Watt.dosage(dos)
#' @param dos  a ni X nl dosage matrix containing the number of derived/alternate alleles each individual carries
#' at each SNP
#' @param L the length of the sequence
#' @return if \code{L=NULL} (default), returns \eqn{\theta_{Watterson}}, else return \eqn{\theta_{Watterson}/L}
#' 
################################################### 
#' @export

theta.Watt.dosage<-function(dos,L=NULL){
if(!is.matrix(dos)){
  if(class(dos)[[1]]=="bed.matrix") dos<-gaston::as.matrix(dos) 
  else dos <- as.matrix(dos)
}
lims <- range(dos, na.rm = TRUE)
if ((lims[2] > 2) | (lims[1] < 0)) 
  stop("input dosage matrix should contains only 0, 1 and 2s")

nis<-nrow(dos)  
p<-colMeans(dos,na.rm=TRUE)/2
SegSites<-sum(p!=0.0 & p!=1.0,na.rm=TRUE)
a<-sum(1/(1:(nis*2-1)))
if(is.null(L)) SegSites/a else 1/L*SegSites/a
}


#####################################################
#' @title Estimates Tajima's D
#' @description Estimates Tajima's D from dosage data
#' @usage TajimaD.dosage(dos)
#' @param dos  a ni X nl dosage matrix containing the number of derived/alternate alleles each individual carries
#' at each SNP
#' @return Tajima's D (eqn 38 of Tajima, 1989)
#' @references
#' 
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203831/}{Tajima F. 1989} 
#' Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism. Genetics 123:585–595. 
#' 
#' @example 
#' 
#' 
#####################################################
#' @export


TajimaD.dosage<-function(dos){
  if(!is.matrix(dos)){
    if(class(dos)[[1]]=="bed.matrix") dos<-gaston::as.matrix(dos) 
    else dos <- as.matrix(dos)
  }
  lims <- range(dos, na.rm = TRUE)
  if ((lims[2] > 2) | (lims[1] < 0)) 
    stop("input dosage matrix should contains only 0, 1 and 2s")

  nis<-2*nrow(dos)  
  
  p<-colMeans(dos,na.rm=TRUE)/2
  Ss<-sum(p>0.0 & p<1.0,na.rm=TRUE)
  a<-sum(1/(1:(nis-1)))
  a2<-sum(1/((1:(nis-1))^2))
  tW<-Ss/a
  pis<-pi.dosage(dos)
  tmp1<- ((nis+1)/(3*(nis-1))-1/a)/a
  num.tmp2<-2*(nis^2+nis+3)/(9*nis*(nis-1))-(nis+2)/(nis*a)+a2/(a^2)
  den.tmp2<-a^2+a2
  VTD<-tmp1*Ss+num.tmp2/den.tmp2*Ss*(Ss-1)
  (pis-tW)/sqrt(VTD)
  
}
  
  
