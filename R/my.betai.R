betai <- function(dat){
  # estimates betas from Weir & Hill 2002 Ann Rev Genet
  # gendata is a data frame with first column containing pop id and remaining cols containing genotypes
  # betai are from (7) of WH2002
  # betaio are population estimates over loci 1-(sum of nums)/(sum of dens)
  # betaw are average over populations using (9) of WH2002
  # betaii' can be extracted with betas[1,,] etc...
  ind.count.n<-function(data){
    dum<-function(x){
      a<-!is.na(x)
      tapply(a,data[,1],sum)
    }
    data[,1]<-factor(data[,1])
    apply(data[,-1],2,dum)
    
  }

  npop<-length(table(dat[,1]))
  nloc<-dim(dat)[2]-1
  al.c.t<-apply(dat[,-1],2,function(x) table(factor(x)))
  
  al.c.pp<-apply(dat[,-1],2,function(x) tapply(factor(x),dat[,1],table))
  al.c.pp.m<-lapply(al.c.pp,function(x) matrix(unlist(x),byrow=TRUE,nrow=npop))
  pb<-lapply(al.c.t,function(x) {dum<-sum(x);if (dum>0.0) x/dum})
  p<-lapply(al.c.pp.m,function(x) {apply(x,1,
                                         function(y) {dum<-sum(y); if (dum>0.0) y/dum else y})
  })
  
  
  ninds <- ind.count.n(dat) 
  n2 <- ninds^2
  n2 <- sweep(n2,2,apply(ninds,2,sum),FUN="/")
  nic <- ninds-n2 #top of p729, WH2002
  snic <- apply(nic,2,sum,na.rm=TRUE) #sum over pop
  betas <- array(numeric(npop^2*nloc),dim=c(npop,npop,nloc))
  nums <- array(numeric(npop^2*nloc),dim=c(npop,npop,nloc))
  lden <- numeric(nloc)
  for (il in 1:nloc){
    dum1 <- sweep(sweep(p[[il]],1,pb[[il]],FUN="-")^2,2,ninds[,il],FUN="*")
    dum2 <- sweep(p[[il]]*(1.0-p[[il]]),2,nic[,il],FUN="*")
    sden <- sum(dum1+dum2,na.rm=TRUE)
    lden[il] <- sden
    for (ip1 in 1:npop){
      for (ip2 in ip1:npop){
        if (ip1==ip2)  
          if (ninds[ip1,il]>0){
            dum1 <- p[[il]][,ip1]*(1.0-p[[il]][,ip1])
            nums[ip1,ip1,il] <- snic[il]*ninds[ip1,il]/(ninds[ip1,il]-1.0)*sum(dum1)
            betas[ip1,ip1,il] <- 1.0-nums[ip1,ip1,il]/sden
          }
        else if ((ninds[ip1,il]>0) & (ninds[ip2,il]>0)){
          dum1 <- p[[il]][,ip1]*(1.0-p[[il]][,ip2])
          dum2 <- p[[il]][,ip2]*(1.0-p[[il]][,ip1])
          nums[ip1,ip2,il] <- snic[il]*sum(dum1+dum2)
          betas[ip1,ip2,il]  <- 1.0-nums[ip1,ip2,il]/2.0/sden
        }
      }
    }
  }
  betai  <- t(apply(betas,3,diag)) # betai per pop and locus
  
  betaio  <- 1-apply(apply(nums,3,diag),1,sum)/sum(lden) # beta per pop over loci, as sums of numerators divided by sums of denominators
  
  betaw  <- 1-apply(apply(nums,3,diag)*ninds,2,sum)/(apply(ninds,2,sum)*lden) # beta overall per loci, eq 9 p 734 of Weir & Hill 2002,
  
  return(list(betaiip=betas,nic=snic,nums=nums,lden=lden,ninds=ninds,nic=nic,betai=betai,betaiov=betaio,betaw=betaw))
}


##############################################################
#' @title Estimation of \eqn{\beta} per population
#' 
#' @description Estimates \eqn{\beta} per population, using unweighted 
#' allele frequencies
#' 
#' @usage betai.unw(dat)
#' 
#' @param dat a dataframe
#' 
#'
#' @return Betai a table (nloc*npop) of individual loci \eqn{\beta_i}
#' @return Betai.ovl a vector of population's \eqn{\beta_i}, average over loci
#' @return BetaW a vector of locus \eqn{\beta_w}, where average is taken for each 
#' locus over populations
#' @return Hi a table (nloc*npop) of within population diversities
#' @return Hb a vector (nloc) of between populations diversities
#' 
#' @details  Beware, only valid for non inbred diploids.
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#'
#' @examples 
#'
#'  dat<-sim.genot(size=100,N=c(100,1000,10000),nbloc=50,nbal=10)
#'  betai.unw(dat)$Betai 
#'
#' @export
#' 
##############################################################
betai.unw<-function(dat){
  ind.count.n<-function(data){
    dum<-function(x){
      a<-!is.na(x)
      tapply(a,data[,1],sum)
    }
    data[,1]<-factor(data[,1])
    apply(data[,-1],2,dum)
  }
  npop<-length(table(dat[,1]))
  nloc<-dim(dat)[2]-1
  al.c.pp<-apply(dat[,-1],2,function(x) tapply(factor(x),dat[,1],table))
  al.c.pp.m<-lapply(al.c.pp,function(x) matrix(unlist(x),byrow=TRUE,nrow=npop))
  
  p<-lapply(al.c.pp.m,function(x) {apply(x,1,
                                         function(y) {dum<-sum(y); if (dum>0.0) y/dum else y})
  })
  
  
  ninds <- ind.count.n(dat) 
  
  H <- array(numeric(npop^2*nloc),dim=c(npop,npop,nloc))
  lden <- numeric(nloc);Hw<-numeric(nloc);Hb<-numeric(nloc)
  acc.p<-numeric(nloc);acc.pp<-numeric(nloc)
  for (il in 1:nloc){
    dum1 <- ninds[,il]/(ninds[,il]-1)*(1-apply(p[[il]]^2,2,sum))
      for (ip1 in 1:npop){
      if (ninds[ip1,il]>1){
        for (ip2 in ip1:npop){
          if(ip2==ip1){
            acc.p[il]<-acc.p[il]+1
            H[ip1,ip1,il]<-dum1[ip1]
            Hw[il]<-Hw[il]+dum1[ip1]
          } else    
        if(ninds[ip2,il]>1){
         acc.pp[il]<-acc.pp[il]+1
         dum2<-1-sum(p[[il]][,ip1]*p[[il]][,ip2])
         H[ip1,ip2,il]<-dum2
         Hb[il]<-Hb[il]+dum2
        }
        }
  }
  }
  }
      Hw<-Hw/acc.p
      Hb<-Hb/acc.pp
      Hi<-apply(H,3,diag)
      Betai<-1-sweep(Hi,2,Hb,FUN="/")
      BetaW<-1-Hw/Hb
      Betai.ovl<-numeric(npop)
      for (ip in 1:npop){
        lt<-which(ninds[ip,]>1)
        Betai.ovl[ip]<-1-sum(Hi[ip,lt])/sum(Hb[lt])
}
    return(list(Betai=Betai,Betai.ovl=Betai.ovl,BetaW=BetaW,Hi=Hi,Hb=Hb))  
  }
