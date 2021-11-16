## ----message=FALSE------------------------------------------------------------
library(hierfstat) #load the library
data(diploid) # info about this data set with ?diploid
head(diploid)

## -----------------------------------------------------------------------------
table(diploid[,1])

## -----------------------------------------------------------------------------
data(cont.isl99)
head(cont.isl99)

## -----------------------------------------------------------------------------
dip<-hierfstat::read.fstat(system.file("extdata","diploid.dat",package="hierfstat"))
head(dip)

## ----message=FALSE------------------------------------------------------------
library(adegenet,quietly=TRUE)
data(nancycats)
head(genind2hierfstat(nancycats)[,1:10]) # only the first 10 loci
data(H3N2)
head(genind2hierfstat(H3N2,pop=rep(1,dim(H3N2@tab)[1]))[,1:10]) # only the first 10 positions
data(eHGDP)
head(genind2hierfstat(eHGDP))[,1:11] 

## ----message=FALSE------------------------------------------------------------
#basic.stats(nancycats)
hierfstat::Hs(nancycats) #mean populations gene diversities
hierfstat::Ho(nancycats) # mean populations observed heterozygosities
#genet.dist(nancycats)

## ----message=FALSE------------------------------------------------------------
obj <-read.snp(system.file("files/exampleSnpDat.snp",package="adegenet"), chunk=2,quiet=TRUE)

as.matrix(obj)


## ----message=FALSE------------------------------------------------------------
library(gaston,quietly=TRUE)
filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
x1 <- read.VCF( filepath )
x1
dim(x1)
with(x1@snps,table(A1,A2))

## ----eval=FALSE---------------------------------------------------------------
#  library(SNPRelate)
#  snpgdsVCF2GDS(filepath, "test1.gds", method="biallelic.only")
#  snpgdsSummary("test1.gds")
#  f<-snpgdsOpen("test1.gds")

## ----eval=FALSE---------------------------------------------------------------
#  x2<-2-snpgdsGetGeno(f)
#  
#  #check that allele frequencies are the same with the two methods
#  all.equal(colMeans(x2)/2,colMeans(as.matrix(x1)/2),check.names=FALSE)

## -----------------------------------------------------------------------------
dat<-qn2.read.fstat(system.file("extdata","qn2_sex.dat",package="hierfstat"))
names(dat)
head(dat$sex)
head(dat$dat[,1:10])
#sexbias.test(dat[[1]],sex=dat[[2]])

## -----------------------------------------------------------------------------
bed<-ms2bed(system.file("extdata","2pops_asspop.txt",package="hierfstat"))
head(as.matrix(bed[,1:10])) #first 10 columns of bed matrix

## -----------------------------------------------------------------------------
fs.dosage(as.matrix(bed),pop=rep(1:2,each=50)) # population specific FSTs

## -----------------------------------------------------------------------------
dat.HGDP<-genind2hierfstat(eHGDP)
dos.HGDP<-fstat2dos(dat.HGDP[,-1])
fs.HGDP<-fs.dosage(dos.HGDP,pop=dat.HGDP[,1])

## -----------------------------------------------------------------------------
eHGDP@other$popInfo[which(fs.HGDP$Fs[2,]>0.15),]
eHGDP@other$popInfo[which(fs.HGDP$Fs[2,]<0.0),]


