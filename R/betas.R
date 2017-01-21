#####################################################
#'
#' @title Estimate \eqn{\beta}s per population and a bootstrap confidence interval
#' 
#' @description Estimate \eqn{\beta}s per population and a bootstrap confidence interval 
#' 
#' @usage betas(dat,nboot=0,lim=c(0.025,0.975),diploid=TRUE,betaijT=FALSE)
#' 
#' @param dat data frame with genetic data and pop identifier
#' @param nboot number of bootstrap samples.
#' @param lim width of the bootstrap confidence interfal
#' @param diploid whether the data comes from a diploid organism
#' @param betaijT whether to estimate individual coancestries
#' 
#' @return Hi Within population gene diversities
#' @return Hb Between populations gene diversities
#' @return betaiovl Average \eqn{\beta_i} over loci
#' @return betaW Average of the betaiovl
#' @return ci The bootstrap confidence interval 
#' (only if more than 100 bootsraps requested and if more than 10 loci are present)
#' @return if betaijT=TRUE, return only the matrix of individual coancestries. 
#' individual inbreeding coeeficient can be obtained by multiplying by 2 the diagonal 
#' and substracting 1 
#' 
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @references \href{http://biorxiv.org/content/early/2016/11/17/088260}{Weir and Goudet } A unified characterization of population structure and relatedness.
#'
#' 
#' @examples
#' \dontrun{
#' #3 different populatoon sizes leads to 3 different betais
#' dat<-sim.genot(size=40,N=c(50,200,1000),nbloc=50,nbal=10)
#' betas(dat,nboot=100)
#'  
#' #individual coancestries from the smallest population are large
#' ind.coan<-betas(cbind(1:120,dat[,-1]),betaij=T)
#' image(1:120,1:120,ind.coan,xlab="Inds",ylab="Inds")
#' 
#' #extracting individual inbreeding coefficients
#' dat<-sim.genot(size=20,nbloc=100,nbal=20,mig=0.01,f=c(0,0.3,0.7))
#' ind.coan<-betas(cbind(1:60,dat[,-1]),betaijT=TRUE)
#' ind.inb<-(diag(ind.coan)*2-1)
#' hist(ind.inb,breaks=seq(-1,1,0.05),xlab="Inidivual inbreeding coeeficients");abline(v=c(0,0.3,0.7),col="red",lwd=2)
#' 
#' }
#'
#'@export
#'
#####################################################
betas<-function(dat,nboot=0,lim=c(0.025,0.975),diploid=TRUE,betaijT=FALSE){
  ratio.Hi.Hb<-function(x){
    dum<-which(!is.na(x))
    sum(x[dum])/sum(Hb[dum]) 
  }
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  
  pfr<-pop.freq(dat,diploid)
  pfr2<-lapply(pfr,function(x) t(x) %*% x)
  nl<-dim(dat)[2]-1
  np<-length(table(dat[,1]))
  ns<-ind.count.n(dat)
  if (diploid) ns<-ns*2
  ns[is.na(ns)]<-0
  Hi<-matrix(numeric(np*nl),ncol=nl)
  dimnames(Hi)<-dimnames(ns)
  Hb<-numeric(nl)
  for (il in 1:nl){
    Hi[,il]<-ns[,il]/(ns[,il]-1)*(1-diag(pfr2[[il]]))
    #   Hi[is.na(Hi)]<-0.0
    npl<-sum(ns[,il]>0,na.rm=TRUE)     
    Hb[il]<-1-1/npl/(npl-1)*sum((pfr2[[il]]-diag(diag(pfr2[[il]]))),na.rm=TRUE)     
  }
  betai<-1-apply(Hi,1,ratio.Hi.Hb)
  betaW<-mean(betai,na.rm=T)
  if (betaijT) {
    ratio.Mij.Hb <- function(x) {
      dum <- which(!is.na(x))
      a<-mean(x[dum])/mean(Hb[dum])
      b<-1/mean(Hb[dum])
      a-b
    }
    H<-array(dim=c(np,np,nl))
    for (il in 1:nl) H[,,il]<-pfr2[[il]]
    betaij<-1+apply(H,c(1,2),ratio.Mij.Hb)
    return(betaij)
  }
  if (nboot==0){return(list(Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW))}
  if (nboot<100){
    warning("Less than 100 bootstrap requested, can't estimate Conf. Int.")
    return(list(Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW))
  }
  else
  {
    if (nl<10) {
      warning("Less than 10 loci, can't estimate Conf. Int.")
      return(list(Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW))
      }
    
    boot.bi<-matrix(numeric(nboot*np),nrow=nboot)
    nls<-apply(ns,1,function(x) which(x>0))
    if (is.matrix(nls)){
    for (ib in 1:nboot){
      for (ip in 1:np){
        dum<-sample(nls[,ip],replace=TRUE)
        boot.bi[ib,ip]<-1-sum(Hi[ip,dum])/sum(Hb[dum])
      }}}
    else{
      for (ib in 1:nboot){
        for (ip in 1:np){
          dum<-sample(nls[[ip]],replace=TRUE)
          boot.bi[ib,ip]<-1-sum(Hi[ip,dum])/sum(Hb[dum])
          
    }}}
    bi.ci<-apply(boot.bi,2,quantile,lim,na.rm=TRUE)
    return(list(Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW,ci=bi.ci))
  }
    
}
