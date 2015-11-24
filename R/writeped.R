#' @title Write ped file plink analysis
#' @description write a ped and a map file suitable for analysis with 
#' \href{plink}{http://pngu.mgh.harvard.edu/~purcell/plink/} 
#' @usage write.ped(dat, ilab = NULL, pop = NULL, 
#'         fname = "dat",na.str="0",f.id=NULL,m.id=NULL,loc.pos=NULL) 
#'                  
#' @param dat 
#' @param ilab
#' @param pop
#' @param fname
#' @param na.str 
#' @param f.id
#' @param m.id
#' @param loc.pos
#'                 
#' @return a map file containing the loci positions 
#' @return a ped file containing genotypes etc... 
#'  
#' @references \href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/}{Purcell etal (2007) PLINK:} A Tool Set for Whole-Genome Association 
#'  and Population-Based Linkage Analyses 81:559-575
#'   
#'    
#'    
#'      
#' @export  
#' 
####################################################################################
write.ped<-function (dat, ilab = NULL, pop = NULL, fname = "dat",na.str="0",f.id=NULL,m.id=NULL,loc.pos=NULL) 
{
    dum <- getal.b(dat[, -1])
    dum[is.na(dum)] <- na.str
    nind <- dim(dum)[1]
    nloc <- dim(dum)[2]
	ddum<-matrix(numeric(nind*nloc*2),nrow=nind)
	al2<-(1:nloc)*2
	al1<-al2-1
	ddum[,al1]<-dum[,,1]
	ddum[,al2]<-dum[,,2]
    popid <- dat[, 1]
    locnames <- paste("L", names(dat)[-1], sep = "")
	mapf<-cbind(0,locnames,0,0)
    write.table(mapf, paste(fname,".map",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	datn<-data.frame(fam.id=popid,ind.id=1:nind,p.id=rep(0,nind),m.id=rep(0,nind),sex=rep(0,nind),pheno=rep(0,nind),ddum)
	write.table(datn,paste(fname,".ped",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
