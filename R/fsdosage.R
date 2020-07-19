#####################################################
#'
#' Estimates F-statistics from dosage data
#' 
#' 
#' Reports individual inbreeding coefficients, Population specific Fsts and Fiss from dosage data 
#' 
#' @aliases fst.dosage fis.dosage
#' @usage fs.dosage(dos, pop, matching = FALSE)
#' 
#' 
#' @param dos either a matrix with snps columns and individuals in rows containing allelic dosage (number [0,1 or 2] of alternate alleles); 
#'     or a square matrix with as many rows and columns as the number of individuals and containing the proportion of matching alleles
#' @param pop a vector containing the identifier of the population to which the individual in the corresponding row belongs
#' @param matching logical:TRUE if dos is a square matrix of allelic matching; FALSE otherwise
#' 
#' @return Fi list of individual inbreeding coefficients, estimated with the reference being the population to which the individual belongs. 
#' @return FsM matrix containing population specific FSTs on the diagonal. The off diagonal elements contains the average of the kinships 
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
#' @export
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
    Fsts <- unlist(lapply(wil, function(x) mean(Mij[x, x], na.rm = TRUE)))
	Mb<-0.0
	mMij<-matrix(numeric(npop^2),ncol=npop)
# 	Mb.min<- 10000
   for (i in 2:npop) {
     p1 <- wil[[i]]
	   for (j in 1:(i-1)){
		    p2 <- wil[[j]]
		    mMij[i,j]<-mMij[j,i]<-mean(Mij[p1,p2],na.rm=TRUE)
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

#' @describeIn fs.dosage Plot function for fs.dosage class
#' @param x a fs.dosage object
#' @param ... further arguments to pass
#' @method plot fs.dosage    
#' @export 
plot.fs.dosage<-function(x,...){
  graphics::par(mfrow=c(2,2))
  graphics::boxplot(x$Fi,xlab="pop",ylab=expression(F[i]),main="Individual Inbreeding coeff.",...)
  np<-dim(x$FsM)[1]
  graphics::image(1:np,1:np,x$FsM,main=expression(F[ST]^{XY}),xlab="Pop. X",ylab="Pop. Y")
  graphics::barplot(x$Fs[1,1:np],xlab="pop",ylab="Fis",main="",...);graphics::abline(h=x$Fs[1,np+1])
  graphics::barplot(x$Fs[2,1:np],xlab="pop",ylab="Fst",main="",...);graphics::abline(h=x$Fs[2,np+1])
  graphics::par(mfrow=c(1,1))
}


#' @describeIn fs.dosage Print function for fs.dosage class 
#' @inheritParams plot.fs.dosage
#' @param digits number of digits to print
#' @method print fs.dosage    
#' @export 
#' 
print.fs.dosage<-function(x,digits=4,...){
  print(round(x$Fs,digits=digits))
  invisible(x)
}


################
#' @rdname fs.dosage 
#' @export
################


 fst.dosage<-function(dos, pop, matching = FALSE){
   fs.dosage(dos,pop,matching)$Fs[2,]
 }


################
#' @rdname fs.dosage 
#' @export
################

 fis.dosage<-function(dos, pop, matching = FALSE){
   fs.dosage(dos,pop,matching)$Fs[1,] 
 }