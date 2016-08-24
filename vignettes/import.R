## ------------------------------------------------------------------------
library(hierfstat) #load the library
data(diploid) # info about this data set with ?diploid
head(diploid)

## ------------------------------------------------------------------------
table(diploid[,1])

## ------------------------------------------------------------------------
data(contisl99)
head(cont.isl99)

## ------------------------------------------------------------------------
dip<-read.fstat(system.file("extdata","diploid.dat",package="hierfstat"))
head(dip)

## ------------------------------------------------------------------------
library(adegenet)
data(nancycats)
head(genind2hierfstat(nancycats)[,1:10]) # only the first 10 loci
#basic.stats(nancycats)
#genet.dist(nancycats)
data(H3N2)
head(genind2hierfstat(H3N2,pop=rep(1,dim(H3N2@tab)[1]))[,1:10]) # only the first 10 positions
#basic.stats(genind2hierfstat(H3N2,pop=rep(1,dim(H3N2@tab)[1])),diploid=FALSE)

## ------------------------------------------------------------------------
dat<-qn2.read.fstat(system.file("extdata","qn2_sex.dat",package="hierfstat"))
names(dat)
head(dat$sex)
head(dat$dat[,1:10])
#sexbias.test(dat[[1]],sex=dat[[2]])

## ----eval=FALSE----------------------------------------------------------
#  msdat<-read.ms.snps()

