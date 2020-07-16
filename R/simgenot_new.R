################################################################################
sim.freq<-function(nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0){
  #allows for different N and f for each population
  #modified param so that it reflects correctly population effective size
  genofreq<-function(freq,fi){
    if (fi==0) {geno.freq<-outer(freq,freq)}
    else {geno.freq<-outer(freq,freq)*(1-fi)+diag(nbal)*diag(freq)*fi}
    return(geno.freq)
  }
  if (nbal>99) stop ("Too many alleles, must be <100. Exiting")
  cl<-match.call()
  pl<-vector("list",nbloc)
  freq<-gtools::rdirichlet(nbloc,rep(1,nbal)) #1=uniform freq has dim nbloc,nbal
  if (length(N)!=nbpop) {
    if (length(N)==1) N<-rep(N,nbpop) 
    else stop("N must be a vector of length nbpop. Exiting.")} 
  if (length(f)!=nbpop){
    if (length(f)==1) f<-rep(f,nbpop) 
    else stop("f must be a vector of length nbpop. Exiting.")}
  param<-outer(4*N/(1+f)*(mig+mut),freq,"*") #verify this [nbpop,nbloc,nbal]
  for (il in 1:nbloc){
    x<-matrix(numeric(nbal*nbpop),nrow=nbpop)
    for (ip in 1:nbpop) x[ip,]<-gtools::rdirichlet(1,param[ip,il,])
    pl[[il]]<-x
  }
  gf<-vector("list",nbloc)
  for (il in 1:nbloc){    
    gf[[il]]<-matrix(numeric(nbal^2*nbpop),ncol=nbpop)
    for (ip in 1:nbpop) gf[[il]][,ip]<-genofreq(pl[[il]][ip,],f[ip])}
  if (nbal<10)
    nfun<-function(x,y) x*10+y
  else
    nfun<-function(x,y) x*100+y
  gn<-as.numeric(outer(1:nbal,1:nbal,nfun))
  return(list(call=cl,fpl=pl,gf=gf,gn=gn))
}
##################
#'@export
###################
#########################################################################
sim.genot<-function(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0){
  a<-sim.freq(nbal,nbloc,nbpop,N,mig,mut,f)
  dat<-data.frame(rep(1:nbpop,each=size),matrix(numeric(nbloc*nbpop*size),ncol=nbloc))
  names(dat)<-c("Pop",paste("loc",1:nbloc,sep="."))
  dumf<-function(x) sample(a$gn,size=size,replace=TRUE,prob=x)
  for (il in 1:nbloc) dat[,il+1]<-as.numeric(apply(a$gf[[il]],2,dumf))
  return(dat)
}
################################################################################
sim.freq.t<-function(nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0,t=100){
  #allows for different N and f for each population
  #modified N->Ne so that it reflects correctly population effective size
  genofreq<-function(freq,fi){
    if (fi==0) {geno.freq<-outer(freq,freq)}
    else {geno.freq<-outer(freq,freq)*(1-fi)+diag(nbal)*diag(freq)*fi}
    return(geno.freq)
  }
  if (nbal>99) stop ("Too many alleles, must be <100. Exiting")
  cl<-match.call()
  pl<-vector("list",nbloc)
  freq<-gtools::rdirichlet(nbloc,rep(1,nbal)) 
  if (length(N)!=nbpop) {
    if (length(N)==1) N<-rep(N,nbpop) 
    else stop("N must be a vector of length nbpop. Exiting.")} 
  if (length(f)!=nbpop){
    if (length(f)==1) f<-rep(f,nbpop) 
    else stop("f must be a vector of length nbpop. Exiting.")}
	xmut<-matrix(rep(1/nbal,nbal*nbpop,nrow=nbpop),nrow=nbpop)

    for (il in 1:nbloc){
	xini<-matrix(rep(freq[il,],nbpop),nrow=nbpop,byrow=TRUE)
    xn<-xini
	Ne<-round(N/(1+f))
    for (it in 1:t) {
	xn<-(xn*(1-mig)+mig*xini)*(1-mut)+mut*xmut
	for (ip in 1:nbpop) xn[ip,]<-rmultinom(1,Ne[ip]*2,xn[ip,])/2/Ne[ip]
}
pl[[il]]<-xn
}

	
  gf<-vector("list",nbloc)
  for (il in 1:nbloc){    
    gf[[il]]<-matrix(numeric(nbal^2*nbpop),ncol=nbpop)
    for (ip in 1:nbpop) gf[[il]][,ip]<-genofreq(pl[[il]][ip,],f[ip])}
  if (nbal<10)
    nfun<-function(x,y) x*10+y
  else
    nfun<-function(x,y) x*100+y
  gn<-as.numeric(outer(1:nbal,1:nbal,nfun))
  return(list(call=cl,fpl=pl,gf=gf,gn=gn))
}
###################################################
################################################################################
sim.freq.FIM.t<-function(nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0,t=100){
  #allows for different N and f for each population
  #modified N->Ne so that it reflects correctly population effective size
  genofreq<-function(freq,fi){
    if (fi==0) {geno.freq<-outer(freq,freq)}
    else {geno.freq<-outer(freq,freq)*(1-fi)+diag(nbal)*diag(freq)*fi}
    return(geno.freq)
  }
  if (nbal>99) stop ("Too many alleles, must be <100. Exiting")
  cl<-match.call()
  pl<-vector("list",nbloc)
  freq<-gtools::rdirichlet(nbloc,rep(1,nbal)) 
  if (length(N)!=nbpop) {
    if (length(N)==1) N<-rep(N,nbpop) 
    else stop("N must be a vector of length nbpop. Exiting.")} 
  if (length(f)!=nbpop){
    if (length(f)==1) f<-rep(f,nbpop) 
    else stop("f must be a vector of length nbpop. Exiting.")}
	xmut<-matrix(rep(1/nbal,nbal*nbpop,nrow=nbpop),nrow=nbpop)

    for (il in 1:nbloc){
	xpool<-matrix(rep(freq[il,],nbpop),nrow=nbpop,byrow=TRUE)
    xn<-xpool
	Ne<-round(N/(1+f))
	tot<-sum(Ne)
    for (it in 1:t) {
	xn<-(xn*(1-mig)+mig*xpool)*(1-mut)+mut*xmut
	for (ip in 1:nbpop) xn[ip,]<-rmultinom(1,Ne[ip]*2,xn[ip,])/2/Ne[ip]
	#xpool<-matrix(rep((N/tot)%*%xn,nbpop),nrow=nbpop,byrow=TRUE) #weighted
	xpool<-matrix(rep((rep(1,nbpop)/nbpop)%*%xn,nbpop),nrow=nbpop,byrow=TRUE) #unweighted
}
pl[[il]]<-xn
}

	
  gf<-vector("list",nbloc)
  for (il in 1:nbloc){    
    gf[[il]]<-matrix(numeric(nbal^2*nbpop),ncol=nbpop)
    for (ip in 1:nbpop) gf[[il]][,ip]<-genofreq(pl[[il]][ip,],f[ip])}
  if (nbal<10)
    nfun<-function(x,y) x*10+y
  else
    nfun<-function(x,y) x*100+y
  gn<-as.numeric(outer(1:nbal,1:nbal,nfun))
  return(list(call=cl,fpl=pl,gf=gf,gn=gn))
}

###################################################
#'
#'  Simulate data from a non equilibrium continent-island model
#'
#'  This function allows to simulate genetic data from a non-equilibrium continent-island
#'  model, where each island can have a different size and a different inbreeding 
#'  coefficient. 
#'  
#'  This function simulates genetic data under the continent-islands model (IIM=TRUE) 
#' or the finite island model (IIM=FALSE).
#'  In the IIM, a continent of 
#'  infinite size sends migrants to islands of finite sizes \eqn{N_i} at a rate 
#' \eqn{m}. Alleles can also mutate to a new state at a rate \eqn{\mu}. Under this model, 
#'  the expected \eqn{F_{STi}, \theta_i}, can be calculated and compared to empirical
#'  estimates.
#'
#' In this model, \eqn{\theta_t} can be written as a function of population size 
#' \eqn{N_i}, migration rate \eqn{m}, mutation rate \eqn{\mu} and \eqn{\theta_{(t-1)}}.  
#' 
#' The rational is as follows:
#'  
#'  With probability \eqn{\frac{1}{N}}, 2 alleles from 2 different individuals in 
#'  the current generation are sampled from the same individual of the previous 
#'  generation:     
#'  
#'  -Half the time, the same allele is drawn from the parent;
#'   
#'  -The other half, two different alleles are drawn, but they are identical in 
#'  proportion \eqn{\theta_{(t-1)}}.
#'   
#'  -With  probability \eqn{1-\frac{1}{N}}, the 2 alleles are drawn from different 
#'  individuals in the previous generation, in which case they are identical in 
#'  proportion \eqn{\theta_{(t-1)}}.
#'  
#'  This holds providing that neither alleles have mutated or migrated.  This is 
#'  the case with probability \eqn{(1-m)^2 \times (1-\mu)^2}.
#'  If an allele is a mutant or a migrant, then its coancestry with another allele 
#'  is 0 in the infinite continent-islands model (it is not the case in the finite island model).
#'
#' Note also that the mutation scheme assumed is the infinite allele (or site) 
#' model.  If the number of alleles is finite (as will be the case in what follows),
#' the corresponding mutation model is the K-allele model and the mutation rate 
#' has to be adjusted to \eqn{\mu'=\frac{K-1}{K}\mu}.
#'
#' Lets substitute \eqn{\alpha} for  \eqn{(1-m)^2 (1-\mu)^2} and \eqn{x} for
#'  \eqn{\frac{1}{2N}}.  
#' 
#' The expectation of \eqn{F_{ST}}, \eqn{\theta} can be written as:
#'
#' \deqn{\theta_t=(\alpha (1-x))^t \theta_0 + \frac{x}{1-x}\sum_{i=1}^t (\alpha (1-x))^i} 
#' 
#'  which reduces to \eqn{\theta_t=\frac{x}{1-x}\sum_{i=1}^t (\alpha (1-x))^i} if \eqn{\theta_0=0}.
#'  
#'  Transition equations for \eqn{theta} in the migrant-pool island model (IIM=FALSE) are given in Rouseet (1996). 
#'  Currently, the migrant pool is made of equal contribution from each island, irrespective of their size.  
#'     
#'
#' @usage sim.genot.t(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,
#' mig=0.001,mut=0.0001,f=0,t=100,IIM=TRUE)
#'  
#' @param size the number of sampled individuals per island
#' @param nbal the number of alleles per locus (maximum of 99)
#' @param nbloc the number of loci to simulate
#' @param nbpop the number of islands to simulate
#' @param N the effective population sizes of each island. If only one number, all
#' islands are assumed to be of the same size
#' @param mig the migration rate from the continent to the islands
#' @param mut the mutation rate of the loci
#' @param f the inbreeding coefficient for each island
#' @param t the number of generation since the islands were created
#' @param IIM whether to simulate a continent island Model (default) or a migrant pool island Model
#' 
#' @return A data frame with size*nbpop rows and nbloc+1 columns. Each row is an
#' individual, the first column contains the island to which the individual belongs, 
#' the following nbloc columns contain the genotype for each locus. 
#'

#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @references 
#'  \href{http://www.ncbi.nlm.nih.gov/pubmed/8846911}{Rousset, F. (1996)} Equilibrium values of measures of population subdivision for 
#'  stepwise mutation processes. Genetics 142:1357
#' 
#' @examples
#' 
#' psize<-c(100,1000,10000,100000,1000000)
#' dat<-sim.genot.t(nbal=4,nbloc=20,nbpop=5,N=psize,mig=0.001,mut=0.0001,t=100)
#' summary(wc(dat)) #Weir and cockerham overall estimators of FST & FIS
#' betas(dat) # Population specific estimator of FST
#' 
#' @export
#' 
#########################################################################
sim.genot.t<-function(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,mig=0.001,mut=0.0001,f=0,t=100,IIM=TRUE){
  if (IIM) a<-sim.freq.t(nbal,nbloc,nbpop,N,mig,mut,f,t)
  else a<-sim.freq.FIM.t(nbal,nbloc,nbpop,N,mig,mut,f,t)
  dat<-data.frame(rep(1:nbpop,each=size),matrix(numeric(nbloc*nbpop*size),ncol=nbloc))
  names(dat)<-c("Pop",paste("loc",1:nbloc,sep="."))
  dumf<-function(x) sample(a$gn,size=size,replace=TRUE,prob=x)
  for (il in 1:nbloc) dat[,il+1]<-as.numeric(apply(a$gf[[il]],2,dumf))
  return(dat)
}
###############################################################################
###############################################################################
###############################################################################
################################################################################
sim.freq.metapop.t<-function(nbal,nbloc,nbpop,N,mig,mut,f,t,epsilon=1e-8){
  #allows for different N and f for each population
  #allows for different migration rates among populations, mig must be a matrix of dim nbpop*nbpop
  #modified param so that it reflects correctly population effective size
  genofreq<-function(freq,fi){
    if (fi==0) {geno.freq<-outer(freq,freq)}
    else {geno.freq<-outer(freq,freq)*(1-fi)+diag(nbal)*diag(freq)*fi}
    return(geno.freq)
  }
  if (nbal>99) stop ("Too many alleles, must be <100. Exiting")
  cl<-match.call()
  pl<-vector("list",nbloc)
  freq<-gtools::rdirichlet(nbloc,rep(1,nbal)) 
  if (length(N)!=nbpop) {
    if (length(N)==1) N<-rep(N,nbpop) 
    else stop("N must be a vector of length nbpop. Exiting.")} 
  if (length(f)!=nbpop){
    if (length(f)==1) f<-rep(f,nbpop) 
    else stop("f must be a vector of length nbpop. Exiting.")}
  #  sanity check for migration matrix
    mig.dim<-dim(mig)
    if (is.null(mig.dim) || (mig.dim[1]!=mig.dim[2]) || (mig.dim[1]!=nbpop) || (abs(sum(rowSums(mig)-1))>epsilon)){
     stop(paste("mig must be a ",nbpop,"X",nbpop," matrix with each row summing to 1. Exiting",sep=""))
      }	
  
  xmut<-matrix(rep(1/nbal,nbal*nbpop,nrow=nbpop),nrow=nbpop)
  
  for (il in 1:nbloc){
    xini<-matrix(rep(freq[il,],nbpop),nrow=nbpop,byrow=TRUE)
    xn<-xini
    tot<-sum(N)
    for (it in 1:t) {
      xn<-(mig%*%xn)*(1-mut)+mut*xmut 
      for (ip in 1:nbpop) xn[ip,]<-rmultinom(1,N[ip]*2,xn[ip,])/2/N[ip]
    }
    pl[[il]]<-xn
  }
  
  
  gf<-vector("list",nbloc)
  for (il in 1:nbloc){    
    gf[[il]]<-matrix(numeric(nbal^2*nbpop),ncol=nbpop)
    for (ip in 1:nbpop) gf[[il]][,ip]<-genofreq(pl[[il]][ip,],f[ip])}
  if (nbal<10)
    nfun<-function(x,y) x*10+y
  else
    nfun<-function(x,y) x*100+y
  gn<-as.numeric(outer(1:nbal,1:nbal,nfun))
  return(list(call=cl,fpl=pl,gf=gf,gn=gn))
}
#########################################################################
###################################################
#'
#'  Simulate genetic data from a metapopulation model
#'
#'  This function allows to simulate genetic data from a metapopulation
#'  model, where each population can have a different size and a different inbreeding 
#'  coefficient, and migration between each population is given in a migration matrix. 
#'  
#'  This function simulates genetic data under a migration matrix model.
#' Each population \eqn{i} sends a proportion of migrant alleles \eqn{m_{ij}} to population \eqn{j} 
#' and receives a proportion of migrant alleles \eqn{m_{ji}} from population \eqn{j}.
#'
#' In this model, \eqn{\theta_t} can be written as a function of population size 
#' \eqn{N_i}, migration rate \eqn{m_{ij}}, mutation rate \eqn{\mu} and \eqn{\theta_{(t-1)}}.  
#' 
#' The rational is as follows:
#'  
#'  With probability \eqn{\frac{1}{N_i}}, 2 alleles from 2 different individuals in 
#'  the current generation are sampled from the same individual of the previous 
#'  generation:     
#'  
#'  -Half the time, the same allele is drawn from the parent;
#'   
#'  -The other half, two different alleles are drawn, but they are identical in 
#'  proportion \eqn{\theta_{(t-1)}}.
#'   
#'  -With  probability \eqn{1-\frac{1}{N_i}}, the 2 alleles are drawn from different 
#'  individuals in the previous generation, in which case they are identical in 
#'  proportion \eqn{\theta_{(t-1)}}.
#'  
#'  This holds providing that neither alleles have mutated or migrated.  This is 
#'  the case with probability \eqn{m_{ii}^2 \times (1-\mu)^2}.
#'  If an allele is a mutant, then its coancestry with another allele 
#'  is 0.
#'
#' Note also that the mutation scheme assumed is the infinite allele (or site) 
#' model.  If the number of alleles is finite (as will be the case in what follows),
#' the corresponding mutation model is the K-allele model and the mutation rate 
#' has to be adjusted to \eqn{\mu'=\frac{K-1}{K}\mu}.
#' 
#' Continue derivation   
#'
#' @usage sim.genot.metapop.t(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,
#' mig=diag(3),mut=0.0001,f=0,t=100)
#'  
#' @param size the number of sampled individuals per population
#' @param nbal the number of alleles per locus (maximum of 99)
#' @param nbloc the number of loci to simulate
#' @param nbpop the number of populations to simulate
#' @param N the effective population sizes of each population. If only one number, all
#' populations are assumed to be of the same size
#' @param mig a matrix with nbpop rows and columns giving the migration rate 
#' from population i (in row) to population j (in column). Each row must sum to 1.
#' @param mut the mutation rate of the loci
#' @param f the inbreeding coefficient for each population
#' @param t the number of generation since the islands were created
#' 
#' @return A data frame with size*nbpop rows and nbloc+1 columns. Each row is an
#' individual, the first column contains the identifier of the population to which the individual belongs, 
#' the following nbloc columns contain the genotype for each locus. 
#'
#'
#'
#' @details   
#'   
#'     
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' @examples
#' 
#' #2 populations
#' psize<-c(10,1000)
#' mig.mat<-matrix(c(0.99,0.01,0.1,0.9),nrow=2,byrow=TRUE)
#' dat<-sim.genot.metapop.t(nbal=10,nbloc=100,nbpop=2,N=psize,mig=mig.mat,mut=0.00001,t=100)
#' betas(dat)$betaiovl # Population specific estimator of FST
#' 
#' #1D stepping stone
#' \dontrun{
#' np<-10
#' m<-0.2
#' mig.mat<-diag(np)*(1-m)
#' diag(mig.mat[-1,-np])<-m/2
#' diag(mig.mat[-np,-1])<-m/2
#' mig.mat[1,1:2]<-c(1-m/2,m/2)
#' mig.mat[np,(np-1):np]<-c(m/2,1-m/2)
#' dat<-sim.genot.metapop.t(nbal=10,nbloc=50,nbpop=np,mig=mig.mat,t=400)
#' pcoa(as.matrix(genet.dist(dat))) # principal coordinates plot
#' }
#' 
#' @export
#' 
##################################################################################################
sim.genot.metapop.t<-function(size=50,nbal=4,nbloc=5,nbpop=3,N=1000,mig=diag(3),mut=0.0001,f=0,t=100){
  a<-sim.freq.metapop.t(nbal,nbloc,nbpop,N,mig,mut,f,t)
  dat<-data.frame(rep(1:nbpop,each=size),matrix(numeric(nbloc*nbpop*size),ncol=nbloc))
  names(dat)<-c("Pop",paste("loc",1:nbloc,sep="."))
  dumf<-function(x) sample(a$gn,size=size,replace=TRUE,prob=x)
  for (il in 1:nbloc) dat[,il+1]<-as.numeric(apply(a$gf[[il]],2,dumf))
  return(dat)
}
