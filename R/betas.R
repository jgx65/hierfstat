#####################################################
#'
#' Estimates \eqn{\beta}s per population and a bootstrap confidence interval
#' 
#' Estimate populations (Population specific FST) or individual coancestries 
#' and a bootstrap confidence interval, assuming random mating 
#' 
#' If betaijT=TRUE, and the first column contains a unique identifier for 
#' each individual, the function returns the matrix of individual coancestries/kinships.  
#' Individual inbreeding coefficients can be obtained by multiplying by 2 the diagonal 
#' and substracting 1.
#' 
#' @usage betas(dat,nboot=0,lim=c(0.025,0.975),diploid=TRUE,betaijT=FALSE)
#' 
#' @param dat data frame with genetic data and pop identifier
#' @param nboot number of bootstrap samples.
#' @param lim width of the bootstrap confidence interval
#' @param diploid whether the data comes from a diploid organism
#' @param betaijT whether to estimate individual coancestries
#' 
#' @return Hi Within population gene diversities (complement to 1 of matching probabilities)
#' @return Hb Between populations gene diversities 
#' @return betaiovl Average \eqn{\beta_i} over loci (Population specific FSTs)
#' @return betaW Average of the betaiovl (overall population FST)
#' @return ci The bootstrap confidence interval of population specific FSTs
#' (only if more than 100 bootstraps requested AND if more than 10 loci are present)
#' @return if betaijT=TRUE, return the matrix of pairwise coancestries only. 
#' 
#' 
#'  
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @references \href{https://www.genetics.org/content/206/4/2085}{Weir and Goudet, 2017 (Genetics)} 
#' A unified characterization of population structure and relatedness.
#'
#' 
#' @examples
#' \dontrun{
#' #3 different population sizes lead to 3 different betais
#' dat<-sim.genot(size=40,N=c(50,200,1000),nbloc=50,nbal=10)
#' betas(dat,nboot=100)
#'  
#' #individual coancestries from the smallest population are large
#' ind.coan<-betas(cbind(1:120,dat[,-1]),betaij=T)
#' graphics::image(1:120,1:120,ind.coan$betaij,xlab="Inds",ylab="Inds")
#' 
#' #extracting individual inbreeding coefficients
#' dat<-sim.genot(size=20,nbloc=100,nbal=20,mig=0.01,f=c(0,0.3,0.7))
#' ind.coan<-betas(cbind(1:60,dat[,-1]),betaijT=TRUE)$betaij
#' ind.inb<-(diag(ind.coan)*2-1)
#' graphics::boxplot(ind.inb~dat[,1], main="Individual inbreeding coefficients")
#' graphics::points(1:3,c(0,0.3,0.7),pch=16,col="red",cex=2)
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
  
  if (betaijT) {inames<-dat[,1];dat[,1]<-1:dim(dat)[1]}
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
 	rownames(betaij)<-colnames(betaij)<-inames
   all.res<-list(betaijT=betaijT,betaij=betaij)
  }
  else {
  if (nboot==0L){
    all.res<-list(betaijT=betaijT,Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW)
  }
  else{
 
  if (nboot<100L){
    warning("Less than 100 bootstrap requested, can't estimate Conf. Int.")
    all.res<-list(betaijT=betaijT,Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW)
  }

    if (nl<10L) {
      warning("Less than 10 loci, can't estimate Conf. Int.")
      all.res<-list(betaijT=betaijT,Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW)
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
    boot.bW<-rowMeans(boot.bi)
    bi.ci<-apply(cbind(boot.bi,boot.bW),2,stats::quantile,lim,na.rm=TRUE)
    
    all.res<-list(betaijT=betaijT,Hi=Hi,Hb=Hb,betaiovl=betai,betaW=betaW,ci=bi.ci)
  }}
  
class(all.res)<-"betas"    
all.res
}

#' @describeIn betas print function for betas class
#' @param x a betas object
#' @param digits number of digits to print
#' @param ... further arguments to pass to print
#' 
#' @method print betas
#' @export
#' 
#' 
####################################################################################
print.betas<-function(x,digits=4,...){
  if (x$betaijT) {cat("   ***Kinship coefficients*** \n \n ")
                      print(round(x$betaij,digits=digits,...))}
  else {
    cat("   ***Population specific and Overall Fst*** \n\n")
    res<-c(x$betaiovl,x$betaW)
    names(res)<-c(names(x$betaiovl),"Overall")
    print(res,digits=digits,...)
    if (!is.null(x$ci)){
      cat("\n   ***Bootstrap CI for population specific and overall FSTs*** \n\n")
      colnames(x$ci)<-names(res)
      print(x$ci,digits=digits,...)
    }
  }
  invisible(x)
}

