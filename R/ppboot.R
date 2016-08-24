################################################################################
pp.sigma.loc<-function(x,y,dat=dat,diploid=TRUE,...){
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  dum<-names(dat)[1]
  wpop<-dat[,1]==x | dat[,1]==y
  ndat<-dat[wpop,]
  ndat[ndat[,1]==x,1]<-1
  ndat[ndat[,1]==y,1]<-2
  return(wc(ndat,diploid,...)$sigma.loc)
}
################################################################################

pp.fst<-function(dat=dat,diploid=TRUE,...){
  cl<-match.call()
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  dat<-dat[order(dat[,1]),]
  Pop<-dat[,1]
  npop<-length(unique(Pop))
  nloc<-dim(dat)[2]-1
  ppsl<-array(numeric(npop*npop*nloc*3),dim=c(npop,npop,nloc,3))
  for (i in 1:(npop-1))
    for (j in (i+1):npop){
      #cat(i," ",j,"\n") #for debugging
      # TODO: optimize by using function for 2 pops only. Needs a new function
      ppsl[i,j,,]<-as.matrix(pp.sigma.loc(unique(Pop)[i],unique(Pop)[j],dat,diploid,...))
    }
  ovl<-apply(ppsl,c(1,2,4),sum)
  ppfst<-ovl[,,1]/apply(ovl,c(1,2),sum)
  res<-list(call=cl,fst.pp=ppfst,vc.per.loc=ppsl)
  class(res)<-"pp.fst"
  res
}

print.pp.fst<-function(x,...){
  print(x$fst.pp)
  invisible(x)
}

################################################################################

boot.ppfst<-function(dat=dat,nboot=100,quant=c(0.025,0.975),diploid=TRUE,...){
  cl<-match.call()
  if (adegenet::is.genind(dat)) dat<-genind2hierfstat(dat)
  typ<-dim(dat)
  if(length(dim(dat))==2){
    #Pop<-dat[,1]
    npop<-length(table(dat[,1]))
    nloc<-dim(dat)[2]-1
    ppsl<-array(numeric(npop*npop*nloc*3),dim=c(npop,npop,nloc,3))
    x<-unique(dat[,1])
    if(is.factor(dat[,1])) dat[,1]<-as.integer(dat[,1])
    for (i in 1:(npop-1)) for (j in (i+1):npop) {
      #cat(i," ",j,"\n") #for debugging
      ppsl[i,j,,]<-as.matrix(pp.sigma.loc(i,j,dat,diploid,...))
    }
  }
  else
  {
    npop<-typ[1]
    nloc<-typ[3]
    ppsl<-dat
  }
  bppfst<-array(numeric(npop*npop*nboot),dim=c(npop,npop,nboot))
  for (i in 1:nboot){
    dum<-sample(nloc,replace=TRUE)
    ovl<-apply(ppsl[,,dum,],c(1,2,4),sum)
    bppfst[,,i]<-ovl[,,1]/apply(ovl,c(1,2),sum)
  }
  #browser()
  ll<-apply(bppfst,c(1,2),quantile,quant[1],na.rm=TRUE)
  hl<-apply(bppfst,c(1,2),quantile,quant[2],na.rm=TRUE)
  dimnames(ll)[[1]]<-dimnames(ll)[[2]]<-x
  dimnames(hl)<-dimnames(ll)
  res<-(list(call=cl,ll=ll,ul=hl,vc.per.loc=ppsl))
  class(res)<-"boot.ppfst"
  res
}

print.boot.ppfst<-function(x,...){
  a<-dim(x$ll)[1]
  ci<-matrix(nrow=a,ncol=a)
  dimnames(ci)<-dimnames(x$ll)
  for (i in 2:a)
    for (j in 1:(a-1)){
      ci[j,i]<-x$ul[j,i]
      ci[i,j]<-x$ll[j,i]
    }
 cat("\n       Upper limit above diagonal \n")
 cat("       Lower limit below diagonal \n \n")
      print(ci)
    invisible(x)
}

################################################################################
boot.ppfis<-function(dat=dat,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4,...){
  cl<-match.call()
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  bs<-basic.stats(dat)
  Ho<-bs$Ho
  Hs<-bs$Hs
  nloc<-dim(Hs)[1]
  npop<-dim(Hs)[2]
  my.boot<-matrix(numeric(nboot*npop),ncol=npop)
  for (i in 1:nboot){
    x<-sample(nloc,replace=TRUE)
    my.boot[i,]<-1-colSums(Ho[x,],na.rm=TRUE)/colSums(Hs[x,],na.rm=TRUE)
  }
  ll<-apply(my.boot,2,quantile,quant[1],na.rm=TRUE)
  hl<-apply(my.boot,2,quantile,quant[2],na.rm=TRUE)
  res<-data.frame(ll=ll,hl=hl)
  return(list(call=cl,fis.ci=round(res,digits=dig)))
}
