#' @title Write ped file for analyses with PLINK
#' @description write a ped and a map file suitable for analysis with 
#' \href{https://www.cog-genomics.org/plink2}{PLINK} 
#' @usage write.ped(dat, ilab = NULL, pop = NULL, 
#'         fname = "dat",na.str="0",f.id=NULL,m.id=NULL,loc.pos=NULL,sex=NULL) 
#' @param dat a hierfstat data frame. if pop=NULL, the first column should contain the population identifier, 
#' otherwise it contains genotypes at the first locus
#' @param ilab individual labels
#' @param pop population id
#' @param fname filename for ped file
#' @param na.str character string to use for missing values
#' @param f.id father id. default to unknown
#' @param m.id mother id. default to unknown
#' @param loc.pos numeric vector of physical positions in base pairs, one per
#' locus in genotype-column order (names are not used for matching). Values
#' must be finite, non-negative whole numbers. NULL writes zero for unknown
#' positions. Positions are written to column four of the MAP file;
#' chromosome and genetic-distance columns remain zero.
#' @param sex the individual sex. default to unknown
#' @return a map file containing the loci positions 
#' @return a ped file containing genotypes etc... 
#' @references \href{https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533}{Chang et al. (2015)} 
#' Second-generation PLINK: rising to the challenge of larger and richer datasets 
#' @export  
####################################################################################
write.ped<-function(dat, ilab = NULL, pop = NULL, fname = "dat",na.str="0",f.id=NULL,m.id=NULL,loc.pos=NULL,sex=NULL) 
{
    if (is.null(pop)) dum <- getal.b(dat[, -1, drop=FALSE]) else dum<-getal.b(dat)
    dum[is.na(dum)] <- na.str
    nind <- dim(dum)[1]
    nloc <- dim(dum)[2]
    if (is.null(loc.pos)) loc.pos <- rep(0, nloc)
    if (!is.numeric(loc.pos) || !is.null(dim(loc.pos)) ||
        length(loc.pos) != nloc || any(!is.finite(loc.pos)) ||
        any(loc.pos < 0) || any(loc.pos != floor(loc.pos)))
        stop("loc.pos must contain one finite non-negative whole number per locus")
	ddum<-matrix(numeric(nind*nloc*2),nrow=nind)
	al2<-(1:nloc)*2
	al1<-al2-1
	ddum[,al1]<-dum[,,1]
	ddum[,al2]<-dum[,,2]
  if (is.null(pop)) popid <- dat[, 1] else popid<-pop
	if (is.null(ilab)) ind.id<-1:nind else ind.id<-ilab
	if(is.null(f.id))  f.id<-rep(0,nind) 
	if(is.null(m.id)) m.id<-rep(0,nind)
	if(is.null(sex)) sex<-rep(0,nind)
	
  if (is.null(pop)) locnames <- paste("L", names(dat)[-1], sep = "") 
	 else locnames <- paste("L", names(dat), sep = "")
	mapf<-cbind(0,locnames,0,format(loc.pos, scientific=FALSE, trim=TRUE))
    utils::write.table(mapf, paste(fname,".map",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	datn<-data.frame(fam.id=popid,ind.id=ind.id,f.id=f.id,m.id=m.id,sex=sex,pheno=rep(0,nind),ddum)
	utils::write.table(datn,paste(fname,".ped",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
