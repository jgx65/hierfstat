################
#' @export
################
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
################
#' @export
################
getal<-function(data){
        #transform a table of genotypes into a table of alleles
		#following anders Goncalves suggestion, replaced nbpop<-max(dat[,1]) with length(unique(dat[,1]))
		#caught exception when encoding alleles with more than 3 digits
		#caught exception when encoding alleles with more than 3 digits but only using allele 1-9
    #added sort to loop following simon forsberg email 
if (is.genind(data)) data<-genind2hierfstat(data)
nc<-dim(data)[2]
true.loci<-2:nc
if (names(data)[nc]=="dummy.loc") true.loci<-2:(nc-1) 
x<-dim(data)
if (max(data[,true.loci],na.rm=TRUE)>1000000) stop("allele encoding with 3 digits maximum")
if (max(data[,true.loci],na.rm=TRUE)<1000000) modulo<-1000
if (max(data[,true.loci],na.rm=TRUE)<10000) {
if (min(data[,true.loci]%/%100,na.rm=TRUE)>=10 & max(data[,true.loci]%%100,na.rm=TRUE)<10) modulo<-1000 else modulo<-100
}
if (max(data[,true.loci],na.rm=TRUE)<100) modulo<-10
firstal<-data[,-1] %/% modulo
secal<-data[,-1] %% modulo
ind<-NULL
nbpop <- length(unique(data[,1]))
for (i in sort(unique(data[,1]))) {
dum<-1:sum(data[,1]==i)
ind<-c(ind,dum)
}
ind<-rep(ind,2)
if (x[2]==2)  data.al<-data.frame(pop=rep(data[,1],2),ind=ind,al=c(firstal,secal))
else data.al<-data.frame(pop=rep(data[,1],2),ind=ind,al=rbind(firstal,secal))
names(data.al)[-c(1:2)]<-names(data)[-1]
return(data.al)
}
#########################################################################
################
#' @export
################
getal.b<-function(data){
  if (is.genind(data)) data<-genind2hierfstat(data)[,-1]
  
x<-dim(data)
if (max(data,na.rm=TRUE)>1000000) stop("allele encoding with 3 digits maximum")
if (max(data,na.rm=TRUE)<1000000) modulo<-1000
if (max(data,na.rm=TRUE)<10000) {
if (min(data%/%100,na.rm=TRUE)>=10 & max(data[,2]%%100,na.rm=TRUE)<10) modulo<-1000 else modulo<-100
}
if (max(data,na.rm=TRUE)<100) modulo<-10
firstal<-data %/% modulo
secal<-data %% modulo
y<-array(dim=c(x,2))
y[,,1]<-as.matrix(firstal)
y[,,2]<-as.matrix(secal)
return(y)
}

#########################################################################
################
#' @export
################
pop.freq<-function(dat,diploid=TRUE){
  lapply(allele.count(dat,diploid),function(x) sweep(x,2,colSums(x,na.rm=TRUE),FUN="/"))
}
#########################################################################
################

################
pop.freq.o<-function(data,diploid=TRUE)
{
if (is.genind(data)) data<-genind2hierfstat(data)
nbloc<-dim(data)[2]-1
nbind<-dim(data)[1]
dumal<-9
nbal<-as.vector(unique(nb.alleles(data.frame(rep(1,nbind),data[,-1]))))
if (length(nbal)==1 & nbal[1]==dumal) dumal<-8
if (diploid) {
data<-data.frame(data,dummy.loc=(sample(dumal,replace=TRUE,size=dim(data)[1])+100)*1001) #to ensure proper output
data<-getal(data)[,-2]
}
else{
data<-data.frame(data,dummy.loc=sample(dumal,replace=TRUE,size=dim(data)[1])+100)
}
freq<-function(x){
#factor(x) necessary for funny allele encoding, but DOES slow down things
  if (is.character(x)) dum<-table(factor(x),data[,1]) else dum<-(table(x,data[,1]))
  eff<-colSums(dum,na.rm=TRUE)
  sweep(dum,2,eff,FUN="/")
}
ndat<-data[,-1]
all.freq<-apply(ndat,2,freq)
if (is.list(all.freq)) all.freq<-all.freq[-(nbloc+1)]
else stop("error in frequency estimation. Exiting")
return(all.freq)
}
#########################################################################
################
#' @export
################
ind.count<-function(data){
 dum<-function(x){
  a<-which(!is.na(x))
  tapply(x[a],data[a,1],length)
 }
 if (is.genind(data)) data<-genind2hierfstat(data)
 
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
    if (is.genind(data)) data<-genind2hierfstat(data)
    
    data[,1]<-factor(data[,1])
    apply(data[,-1],2,dum)
    
  }

#########################################################################
################
#' @export
################
allele.count<-function(data,diploid=TRUE)
{
  if (is.genind(data)) data<-genind2hierfstat(data)
  
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
################
#' @export
################
nb.alleles<-function(data,diploid=TRUE){
if (is.genind(data)) data<-genind2hierfstat(data)
  
dum<-allele.count(data,diploid)
if(is.list(dum))
res<-lapply(dum,fun<-function(x) apply(x,2,fun2<-function(x) sum(x>0)))
else res<-apply(dum,2,fun2<-function(x) sum(x>0))
res<-matrix(unlist(res),nrow=dim(data)[2]-1,byrow=TRUE)
rownames(res)<-names(data)[-1]
return(res)
}
########################################################################
################
#' @export
################
allelic.richness<-function (data, min.n = NULL, diploid = TRUE) 
{
    raref <- function(x) {
        nn <- sum(x)
        dum <- choose(nn - x, min.n)/choose(nn, min.n)
        dum[is.na(dum)] <- 0
        return(sum(1 - dum))
    }
    if (is.genind(data)) data<-genind2hierfstat(data)
    if (dim(data)[2] == 2) 
      data <- data.frame(data, dummy.loc = data[, 2])
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
    Ar<-data.frame(Ar)
    rownames(Ar) <- names(data)[-1]
    colnames(Ar)<-names(table(data[,1]))
	mynas<-which(is.na(t(ind.count(data))),arr.ind=TRUE)
	Ar[mynas]<-NA
    return(list(min.all = min.n, Ar = Ar))
}
#########################################################################
################
#' @export
################
basic.stats<-function(data,diploid=TRUE,digits=4){
# TODO : define plot functions for basic.stats
if (is.genind(data)) data<-genind2hierfstat(data)
if (length(table(data[,1]))<2) data[dim(data)[1]+1,1]<-data[dim(data)[1],1]+1
if(dim(data)[2]==2) data<-data.frame(data,dummy.loc=data[,2])
loc.names<-names(data)[-1]
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
if(!diploid){
#  res[,-2]<-NA
  overall[c(1,9)]<-NA
}
all.res<-list(n.ind.samp=n,pop.freq=lapply(p,round,digits),
Ho=round(sHo,digits),Hs=round(Hs,digits),Fis=round(Fis,digits),
perloc=round(res,digits),overall=round(overall,digits))
class(all.res)<-"basic.stats"
all.res
}

#' @method print basic.stats
#' @export 
print.basic.stats<-function(x,...){
print(list(perloc=x$perloc,overall=x$overall))
invisible(x)
}

################
#' @rdname basic.stats 
#' @export
################
Hs<-function(data,...){
  colMeans(basic.stats(data,...)$Hs,na.rm=TRUE)
}

################
#' @rdname basic.stats 
#' @export
################

Ho<-function(data,...){
  colMeans(basic.stats(data,...)$Ho,na.rm=TRUE)
}  

#####################################################


is.genind<-function (x) 
{
  res <- (methods::is(x, "genind") & methods::validObject(x))
  return(res)
}
