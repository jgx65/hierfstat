#####################################################
#'
#' Estimates F-statistics from dosage data
#' 
#' 
#' Reports individual inbreeding coefficients, Population specific Fsts and Fiss from dosage data 
#' 
#' @usage fs.dosage(dos, pop, matching = FALSE)
#' 
#' 
#' @param dos either a matrix with snps columns and individuals in rows containing allelic dosage (number [0,1 or 2] of alternate alleles); 
#'     or a square matrix with as many rows and columns as the number of individuals and containing the proportion of matching alleles
#' @param pop a vector containing the identifier of the population to which the individual in the corresponding row belongs
#' @param matching logical:TRUE if dos is a square matrix of allelic matching; FALSE otherwise
#' 
#' @return Fi list of individual inbreeding coefficients, estimated with the reference being the population to which the individual belongs. 
#' @return Fsm matrix containing population specific FSTs on the diagonal. The off diagonal elements contains the average of the kinships 
#' for pairs of individuals, one from each population, relative to the mean kinship for pairs of individuals between populations. 
#' @return Fs The first row contains population specific and overall Fis, the second row population specific and overall Fst.
#' 
#'  
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' 
#' \href{http://www.genetics.org/content/206/4/2085}{Weir, BS and Goudet J. 2017} A Unified Characterization 
#' of Population Structure and Relatedness. Genetics (2017) 206:2085 
#'
#' 
#' @examples
#' \dontrun{
#'  dos<-matrix(sample(0:2,size=10000,replace=TRUE),ncol=100)
#'  fs.dosage(dos,pop=rep(1:5,each=20))
#'  plot(fs.dosage(dos,pop=rep(1:5,each=20)))
#' }
#'
#'@export
#######################################
fs.dosage<-function (dos, pop, matching = FALSE)
 #estimates Fis and Fst per population and overall
 
{
    if (!matching) {
        tmp<-beta.dosage(dos,inb=FALSE,Mb=TRUE)
        Mij<-with(tmp,betas*(1-MB)+MB)
        Mii <- (diag(Mij)) * 2 - 1
        diag(Mij) <- NA
    }
    else {
        Mij <- dos
        Mii <- diag(Mij)*2-1
        diag(Mij) <- NA
    }
    pop <- factor(pop)
    x <- levels(pop)
	npop<-length(x)
    wil <- lapply(x, function(z) which(pop == z))
	Fi<-lapply(wil,function(x) Mii[x])
    Fsts <- unlist(lapply(wil, function(x) mean(Mij[x, x], na.rm = T)))
	Mb<-0.0
	mMij<-matrix(numeric(npop^2),ncol=npop)
# 	Mb.min<- 10000
   for (i in 2:npop) {
     p1 <- wil[[i]]
	   for (j in 1:(i-1)){
		    p2 <- wil[[j]]
		    mMij[i,j]<-mMij[j,i]<-mean(Mij[p1,p2])
		    Mb<-Mb+mMij[i,j]
		    #		Mb.min<-min(Mb.min,mean(Mij[p1,p2]))
    }}
	
    Mb <- Mb*2/(npop*(npop-1)) #mean(Mij, na.rm = T)
    diag(mMij)<-Fsts
    resM<-(mMij-Mb)/(1-Mb)
    rownames(resM)<-colnames(resM) <- levels(pop)
    #Fsts <- (Fsts - Mb)/(1 - Mb)
	for (i in 1:npop) Fi[[i]]<-(Fi[[i]] -Fsts[[i]])/(1-Fsts[[i]])
	names(Fi)<-as.character(x)
	Fis<-unlist(lapply(Fi,mean,na.rm=TRUE))
	Fis<-c(Fis,mean(Fis,na.rm=TRUE))
	
    Fsts<- c((Fsts-Mb)/(1-Mb), mean((Fsts-Mb)/(1-Mb), na.rm = TRUE)) #res.Mb
#	res.Mb.min<-c((Fsts-Mb.min)/(1-Mb.min),mean((Fsts-Mb.min)/(1-Mb.min),na.rm=TRUE))
#	res<-data.frame(rbind(res.Mb,res.Mb.min))
    res<-rbind(Fis,Fsts) 
    colnames(res) <- c(levels(pop), "All")
	rownames(res)<-c("Fis","Fst")
    all.res<-(list(Fi=Fi,FsM=resM,Fs=res))
    class(all.res)<-"fs.dosage"
    all.res
}

#' Print function for fs.dosage 
#' 
#' @usage print.fs.dosage(x,...)
#'  
#' @param x a fs.dosage object
#' @param ... further arguments to pass 
#'    
#' @method print fs.dosage    
#' @export 
#' 
print.fs.dosage<-function(x,...){
  print(x$Fs)
  invisible(x)
}

#' Plot function for fs.dosage 
#' 
#' @usage plot.fs.dosage(x,...)
#'  
#' @param x a fs.dosage object
#' @param ... further arguments to pass 
#'    
#' @method plot fs.dosage    
#' @export 
#' 
plot.fs.dosage<-function(x,...){
  par(mfrow=c(2,2))
  boxplot(x$Fi,xlab="pop",ylab=expression(F[i]),main="Individual Inbreeding coeff.")
  np<-dim(x$FsM)[1]
  image(1:np,1:np,x$FsM,main=expression(F[ST]^{XY}),xlab="Pop. X",ylab="Pop. Y")
  barplot(x$Fs[1,],xlab="pop",ylab="Fis",main="")
  barplot(x$Fs[2,],xlab="pop",ylab="Fst",main="")
  par(mfrow=c(1,1))
}


################
#' Reports population Fsts from dosage data
#' 
#' Reports population Fsts from dosage data. see \link{fs.dosage}
#' 
#' @usage fst.dosage(dos, pop, matching = FALSE)
#' 
#' @param dos either a matrix with snps columns and individuals in rows containing allelic dosage (number [0,1 or 2] of alternate alleles); 
#'     or a square matrix with as many rows and columns as the number of individuals and containing the proportion of matching alleles
#' @param pop a vector containing the identifier of the population to which the individual in the corresponding row belongs
#' @param matching logical:TRUE if dos is a square matrix of allelic matching; FALSE otherwise
#' 
#' @return  population specific and overall Fst.
#' @export
################


 fst.dosage<-function (dos, pop, matching = FALSE){
   fs.dosage(dos,pop,matching)$Fs[2,]
 }


################
#'  Reports population Fiss from dosage data
#'
#'Reports population Fiss from dosage data. see \link{fs.dosage}
#'
#'  @usage fis.dosage(dos, pop, matching = FALSE)
#' 
#' @param dos either a matrix with snps columns and individuals in rows containing allelic dosage (number [0,1 or 2] of alternate alleles); 
#'     or a square matrix with as many rows and columns as the number of individuals and containing the proportion of matching alleles
#' @param pop a vector containing the identifier of the population to which the individual in the corresponding row belongs
#' @param matching logical:TRUE if dos is a square matrix of allelic matching; FALSE otherwise
#' 
#' @return  population specific and overall Fis.

#' @export
################

 fis.dosage<-function(dos, pop, matching = FALSE){
   fs.dosage(dos,pop,matching)$Fs[1,] 
 }