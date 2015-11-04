"read.fstat" <-
function (fname, na.s = c("0","00","000","0000","00000","000000","NA"))
{
x<-scan(fname,n=4)
nloc<-x[2]
lnames<-scan(fname,what=character(),skip=1,nlines=nloc)
lnames<-c("Pop",lnames)
dat<-scan(fname,skip=nloc+1,na.strings=na.s)
dat<-data.frame(matrix(dat,ncol=nloc+1,byrow=TRUE))
names(dat)<-lnames
return(dat)
}
############################################################################################
getal<-function(data){
        #transform a table of genotypes into a table of alleles
		#following anders Goncalves suggestion, replaced nbpop<-max(dat[,1]) with length(unique(dat[,1]))
		#caught exception when encoding alleles with more than 3 digits
		#caught exception when encoding alleles with more than 3 digits but only using allele 1-9
    #added sort to loop following simon forsberg email 
x<-dim(data)
if (max(data[,2],na.rm=TRUE)>1000000) stop("allele encoding with 3 digits maximum")
if (max(data[,2],na.rm=TRUE)<1000000) modulo<-1000
if (max(data[,2],na.rm=TRUE)<10000) {
if (min(data[,2]%/%100,na.rm=TRUE)>=10 & max(data[,2]%%100,na.rm=TRUE)<10) modulo<-1000 else modulo<-100
}
if (max(data[,2],na.rm=TRUE)<100) modulo<-10
firstal<-data[,-1] %/% modulo
secal<-data[,-1] %% modulo
ind<-vector(length=0)
nbpop <- length(unique(data[,1]))
for (i in sort(unique(data[,1]))) {
dum<-1:sum(data[,1]==i)
if (i==1) ind<-dum else ind<-c(ind,dum)
}
ind<-rep(ind,2)
if (x[2]==2)  data.al<-data.frame(pop=rep(data[,1],2),ind=ind,al=c(firstal,secal))
else data.al<-data.frame(pop=rep(data[,1],2),ind=ind,al=rbind(firstal,secal))
names(data.al)[-c(1:2)]<-names(data)[-1]
return(data.al)
}
#########################################################################
getal.b<-function(data){
x<-dim(data)
if (max(data[,2],na.rm=TRUE)>1000000) stop("allele encoding with 3 digits maximum")
if (max(data[,2],na.rm=TRUE)<1000000) modulo<-1000
if (max(data[,2],na.rm=TRUE)<10000) {
if (min(data[,2]%/%100,na.rm=TRUE)>=10 & max(data[,2]%%100,na.rm=TRUE)<10) modulo<-1000 else modulo<-100
}
if (max(data[,2],na.rm=TRUE)<100) modulo<-10
firstal<-data %/% modulo
secal<-data %% modulo
y<-array(dim=c(x,2))
y[,,1]<-as.matrix(firstal)
y[,,2]<-as.matrix(secal)
return(y)
}
#########################################################################
pop.freq<-function(data,diploid=TRUE)
{
nbpop<-length(unique(data[,1]))
nbloc<-dim(data)[2]-1
if (diploid) {
data<-data.frame(data,dummy.loc=rep(101,dim(data)[1])*1001) #to ensure proper output
data<-getal(data)[,-2]
}
else{
data<-data.frame(data,dummy.loc=rep(101,dim(data)[1]))
}
freq<-function(x){
#factor(x) necessary for funny allele encoding, but DOES slow down things
  if (is.character(x)) dum<-table(factor(x),data[,1]) else dum<-(table(x,data[,1]))
  eff<-colSums(dum,na.rm=TRUE)
  freq<-sweep(dum,2,eff,FUN="/")
}
ndat<-data[,-1]
all.freq<-apply(ndat,2,freq)
all.freq<-all.freq[-(nbloc+1)]
return(all.freq)
}
#########################################################################
ind.count<-function(data){
 dum<-function(x){
  a<-which(!is.na(x))
  tapply(x[a],data[a,1],length)
 }
 data[,1]<-factor(data[,1])
 apply(data[,-1],2,dum)
}
#########################################################################
  ind.count.n<-function(data){
  #should replace ind.count for ill behaved loci, e.g. many weird alleles
    dum<-function(x){
      a<-!is.na(x)
      tapply(a,data[,1],sum)
    }
    data[,1]<-factor(data[,1])
    apply(data[,-1],2,dum)
    
  }

#########################################################################
allele.count<-function(data,diploid=TRUE)
{
nbpop<-length(unique(data[,1]))
if (diploid) data<-getal(data)[,-2]
nbloc<-dim(data)[2]-1
fun<-function(x){
  dum<-table(x,data[,1])
}
all.count<-apply(as.data.frame(data[,-1]),2,fun)
if (!is.list(all.count))
{
my.dim<-dim(table(data[,2],data[,1]))
my.names<-dimnames(all.count)[[2]]
dum<-list(nbloc)
for (i in seq_along(my.names)){
dum[[i]]<-matrix(all.count[,i],nrow=my.dim[1])
}
names(dum)<-my.names
all.count<-dum
}

return(all.count)
}
#########################################################################
nb.alleles<-function(data,diploid=TRUE){
dum<-allele.count(data,diploid)
if(is.list(dum))
res<-lapply(dum,fun<-function(x) apply(x,2,fun2<-function(x) sum(x>0)))
else res<-apply(dum,2,fun2<-function(x) sum(x>0))
res<-matrix(unlist(res),nrow=dim(data)[2]-1,byrow=TRUE)
rownames(res)<-names(data)[-1]
return(res)
}
########################################################################
allelic.richness<-function (data, min.n = NULL, diploid = TRUE) 
{
    raref <- function(x) {
        nn <- sum(x)
        dum <- choose(nn - x, min.n)/choose(nn, min.n)
        dum[is.na(dum)] <- 0
        return(sum(1 - dum))
    }
    nloc <- dim(data[, -1])[2]
    all.count <- allele.count(data, diploid)
    if (is.null(min.n)) {
        min.n <- 2 * min(ind.count(data), na.rm = TRUE)
        if (!diploid) {
            min.n <- min.n/2
        }
    }
    a <- lapply(all.count, fun1 <- function(x) apply(x, 2, raref))
    Ar <- matrix(unlist(a), nrow = nloc, byrow = TRUE)
    rownames(Ar) <- names(data)[-1]
	mynas<-which(is.na(t(ind.count(data))))
	Ar[mynas]<-NA
    return(list(min.all = min.n, Ar = Ar))
}
#########################################################################
basic.stats<-function(data,diploid=TRUE,digits=4){
# TODO : define plot functions for basic.stats
loc.names<-names(data)[-1]

if (length(table(data[,1]))<2) data[dim(data)[1]+1,1]<-data[dim(data)[1],1]+1
if(dim(data)[2]==2) data<-data.frame(data,dummy.loc=data[,2])
p<-pop.freq(data,diploid)
n<-t(ind.count(data))
if (diploid){
dum<-getal.b(data[,-1])
Ho<-dum[,,1]==dum[,,2]
#sHo<-(n-t(apply(Ho,2,fun<-function(x) tapply(x,data[,1],sum,na.rm=TRUE))))/n
sHo<-(1-t(apply(Ho,2,fun<-function(x) tapply(x,data[,1],mean,na.rm=TRUE))))

mHo<-apply(sHo,1,mean,na.rm=TRUE)
}
else {
sHo<-NA
mHo<-NA
}
sp2<-lapply(p,fun<-function(x) apply(x,2,fun2<-function(x) sum(x^2)))
sp2<-matrix(unlist(sp2),nrow=dim(data[,-1])[2],byrow=TRUE)
if (diploid) {
Hs<-(1-sp2-sHo/2/n)
Hs<-n/(n-1)*Hs
Fis=1-sHo/Hs
}
else {
Hs<-n/(n-1)*(1-sp2)
Fis<-NA
}
np<-apply(n,1,fun<-function(x) sum(!is.na(x)))
mn<-apply(n,1,fun<-function(x) {np<-sum(!is.na(x)); np/sum(1/x[!is.na(x)])}) 
msp2<-apply(sp2,1,mean,na.rm=TRUE)
mp<-lapply(p,fun<-function(x) apply(x,1,mean,na.rm=TRUE))
mp2<-unlist(lapply(mp,fun1<-function(x) sum(x^2)))
if (diploid){
mHs<-mn/(mn-1)*(1-msp2-mHo/2/mn)
Ht<-1-mp2+mHs/mn/np-mHo/2/mn/np
mFis=1-mHo/mHs
}
else{
mHs<-mn/(mn-1)*(1-msp2)
Ht<-1-mp2+mHs/mn/np
mFis<-NA
}
Dst<-Ht-mHs
Dstp<-np/(np-1)*Dst
Htp=mHs+Dstp
Fst=Dst/Ht
Fstp=Dstp/Htp
Dest<-Dstp/(1-mHs)#*np/(np-1)
res<-data.frame(cbind(mHo,mHs,Ht,Dst,Htp,Dstp,Fst,Fstp,mFis,Dest))
names(res)<-c("Ho","Hs","Ht","Dst","Htp","Dstp","Fst","Fstp","Fis","Dest")
if (diploid){
rownames(sHo)<-loc.names
rownames(Fis)<-loc.names
}
is.na(res)<-do.call(cbind,lapply(res,is.infinite))
overall<-apply(res,2,mean,na.rm=TRUE)
overall[7]<-overall[4]/overall[3]
overall[8]<-overall[6]/overall[5]
overall[9]<-1-overall[1]/overall[2]
overall[10]<-overall[6]/(1-overall[2])
names(overall)<-names(res)
all.res<-list(n.ind.samp=n,pop.freq=lapply(p,round,digits),
Ho=round(sHo,digits),Hs=round(Hs,digits),Fis=round(Fis,digits),
perloc=round(res,digits),overall=round(overall,digits))
class(all.res)<-"bas.stats"
all.res

}
print.bas.stats<-function(x,...){
print(list(perloc=x$perloc,overall=x$overall))
invisible(x)
}
#########################################################################
wc<-function(ndat,diploid=TRUE,pol=0.0,trim=FALSE){

#trim is to output only sigloc
  
#added argument pol to specify level of polymorhism wanted
#if loci less polymorphic than pol, then they are left out of the calculations
#setting pol=0.0 will use all loci
# the pol argument is not yet used
#need to modify function to properly handle haploid data
cl<-match.call()
if (!diploid) {
dum<-ndat[,-1]
nd<-max(dum,na.rm=T)
modu<-1000
if (nd<10) modu<-10
if(nd<100) modu<-100
dum<-dum*modu+dum
ndat<-data.frame(ndat[,1],dum)
}

#bs<-basic.stats(ndat,diploid)
pop<-as.integer(ndat[,1])
ni<-length(pop)
#polym<-which(bs$perloc$Ht>=pol,arr.ind=TRUE)
dat<-ndat#[,c(1,polym+1)]
loc.names<-names(dat)[-1]
n<-t(ind.count(dat))
nt<-apply(n,1,sum,na.rm=TRUE)
untyped.loc<-which(nt==0)
typed.loc<-which(nt!=0)

if (length(untyped.loc)>0){
  dat<-dat[,-(untyped.loc+1)]
  n<-t(ind.count(dat))
  nt<-apply(n,1,sum,na.rm=TRUE)
}

alploc<-nb.alleles(cbind(rep(1,ni),dat[,-1]))
np<-dim(n)[2]
npl<-apply(n,1,tempfun<-function(x) sum(!is.na(x))) #accounts for pops not typed for one locus
nl<-dim(n)[1]
p<-pop.freq(dat,diploid)
pb<-pop.freq(cbind(rep(1,length(pop)),dat[,-1]),diploid)
n<-matrix(unlist(n),ncol=np)
nal<-n[rep(1:nl,alploc),]
nc<-(nt-apply(n^2,1,sum,na.rm=TRUE)/nt)/(npl-1)
ntal<-rep(nt,alploc)
ncal<-rep(nc,alploc)
p<-matrix(unlist(lapply(p,t)),ncol=np,byrow=TRUE)
pb<-matrix(unlist(pb),ncol=1)


if (diploid){
dum<-getal.b(dat[,-1])
all.loc<-apply(dum,2,tempfun1<-function(y) as.numeric(dimnames(table(y))[[1]]))

hetpl<-apply(dum,2,fun<-function(z){
lapply(as.numeric(dimnames(table(z))[[1]]),
who.is.het<-function(y) apply(z==y,1,
ind.is.het<-function(x) xor(x[1],x[2])))}
)

mho<-lapply(hetpl,
tempfun2<-function(x) matrix(unlist(lapply(x,
tempfun3<-function(y) tapply(y,pop,sum,na.rm=TRUE))),ncol=np)
)
mho<-matrix(unlist(mho),ncol=np,byrow=TRUE)
mhom<-(2*nal*p-mho)/2
}
else mhom<-nal*p

SSG<-apply(nal*p-mhom,1,sum,na.rm=TRUE)

dum<-nal*(p-2*p^2)+mhom
SSi<-apply(dum,1,sum,na.rm=TRUE)

dum1<-nal*(sweep(p,1,pb))^2
SSP<-2*apply(dum1,1,sum,na.rm=TRUE)

ntalb<-rep(npl,alploc)   #corrects for untyped samples at a locus

MSG<-SSG/ntal
MSP<-SSP/(ntalb-1)
MSI<-SSi/(ntal-ntalb)

sigw<-MSG
sigb<-0.5*(MSI-MSG)
siga<-1/2/ncal*(MSP-MSI)

Fxy<-function(x) x[1]/sum(x,na.rm=TRUE)

FST.pal<-apply(cbind(siga,sigb,sigw),1,Fxy)
FIS.pal<-apply(cbind(sigb,sigw),1,Fxy)

#need to optimize sigloc

loc<-rep(1:nl,alploc)
lsiga<-tapply(siga,loc,sum,na.rm=TRUE)
lsigb<-tapply(sigb,loc,sum,na.rm=TRUE)
lsigw<-tapply(sigw,loc,sum,na.rm=TRUE)

sigloc<-cbind(lsiga,lsigb,lsigw)


if (length(untyped.loc)>0){
x<-order(c(typed.loc,untyped.loc))
dum1<-matrix(numeric(3*length(untyped.loc)),ncol=3)
dum<-rbind(sigloc,dum1,deparse.level=0)[x,]

sigloc<-dum
}
sigloc<-data.frame(sigloc)
names(sigloc)=c("lsiga","lsigb","lsigw")
rownames(sigloc)<-loc.names

lFST<-apply(cbind(lsiga,lsigb,lsigw),1,Fxy)
lFIS<-apply(cbind(lsigb,lsigw),1,Fxy)

tsiga<-sum(siga,na.rm=TRUE)
tsigb<-sum(sigb,na.rm=TRUE)
tsigw<-sum(sigw,na.rm=TRUE)
tFST<-Fxy(c(tsiga,tsigb,tsigw))
tFIS<-Fxy(c(tsigb,tsigw))

res<-list(call=cl,sigma=cbind(loc,siga,sigb,sigw),
sigma.loc=sigloc,
per.al=list(FST=FST.pal,FIS=FIS.pal),
per.loc=list(FST=lFST,FIS=lFIS),
FST=tFST,FIS=tFIS)

if (!diploid){
res<-list(call=cl,sigma=cbind(loc,siga,sigb,sigw),
sigma.loc=sigloc,
per.al=list(FST=FST.pal),
per.loc=list(FST=lFST),
FST=tFST,FIS=NA)
}
class(res)<-"wc"
res
}
print.wc<-function(x,...){
print(list(perloc=x$per.loc,FST=x$FST,FIS=x$FIS))
invisible(x)
}
################################################################################
pp.sigma.loc<-function(x,y,dat=dat,diploid=TRUE,...){
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

boot.ppfst<-function(dat=dat,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4,...){
cl<-match.call()
typ<-dim(dat)
if(length(dim(dat))==2){
#Pop<-dat[,1]
npop<-length(table(dat[,1]))
nloc<-dim(dat)[2]-1
ppsl<-array(numeric(npop*npop*nloc*3),dim=c(npop,npop,nloc,3))
x<-unique(dat[,1])
for (i in x[-npop]) for (j in x[i+1]:x[npop]) {
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
return(list(call=cl,ll=round(ll,digits=dig),ul=round(hl,digits=dig),vc.per.loc=ppsl))
}
################################################################################
boot.ppfis<-function(dat=dat,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4,...){
cl<-match.call()
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
