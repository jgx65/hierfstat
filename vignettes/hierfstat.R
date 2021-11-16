## ----message=FALSE------------------------------------------------------------
library(adegenet)
library(hierfstat)

## -----------------------------------------------------------------------------
data(nancycats) 
is.genind(nancycats)

## -----------------------------------------------------------------------------
bs.nc<-basic.stats(nancycats)
bs.nc
boxplot(bs.nc$perloc[,1:3]) # boxplot of Ho, Hs, Ht

## -----------------------------------------------------------------------------
boxplot(t(allelic.richness(nancycats)$Ar[1:5,])) #5 first loci

## -----------------------------------------------------------------------------
data(gtrunchier) 
wc(gtrunchier[,-2])
varcomp.glob(data.frame(gtrunchier[,1]),gtrunchier[,-c(1:2)])$F #same

## -----------------------------------------------------------------------------
boot.vc(gtrunchier[,1],gtrunchier[,-c(1:2)])$ci

## -----------------------------------------------------------------------------
(Ds<-genet.dist(gtrunchier[,-2],method="Ds")) # Nei's standard genetic distances

## ----eval=FALSE---------------------------------------------------------------
#  pcoa(as.matrix(Ds))

## -----------------------------------------------------------------------------
barplot(betas(nancycats)$betaiovl)


## -----------------------------------------------------------------------------
gt.dos<-fstat2dos(gtrunchier[,-c(1:2)]) #converts fstat format to dosage 
fs.gt<-fs.dosage(gt.dos,pop=gtrunchier[,2]) 
image(1:29,1:29,fs.gt$FsM,main=expression(F[ST]^{XY}))

## -----------------------------------------------------------------------------
x<-indpca(gtrunchier[,-2],ind.labels=gtrunchier[,2])
plot(x,col=gtrunchier[,1],cex=0.5)

## -----------------------------------------------------------------------------
image(beta.dosage(gt.dos),main="kinship and inbreeding \n in Galba truncatula",xlab="",ylab="")

## -----------------------------------------------------------------------------
data("crocrussula")
aic<-AIc(crocrussula$genot)
boxplot(aic~crocrussula$sex)
tapply(aic,crocrussula$sex,mean)
sexbias.test(crocrussula$genot,crocrussula$sex)

