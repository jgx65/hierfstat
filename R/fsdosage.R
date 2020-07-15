 fs.dosage<-function (dos, pop, matching = FALSE)
 #estimates Fis and Fst per population and overall
 
{
    if (!matching) {
        nl <- dim(dos)[2]
        Mij <- 1/2 + tcrossprod(dos - 1)/2/nl
        Mii <- (diag(Mij)) * 2 - 1
        diag(Mij) <- NA
    }
    else {
        Mij <- dos
        Mii <- diag(Mij)
        diag(Mij) <- NA
    }
    pop <- factor(pop)
    x <- levels(pop)
	npop<-length(x)
    wil <- lapply(x, function(z) which(pop == z))
	Fi<-lapply(wil,function(x) Mii[x])
    Fsts <- lapply(wil, function(x) mean(Mij[x, x], na.rm = T))
	Mb<-0.0
# 	Mb.min<- 10000
   for (i in 2:npop) {
	for (j in 1:(i-1)){
        p1 <- wil[[i]]
		p2 <- wil[[j]]
        Mb<-Mb+mean(Mij[p1, p2])
#		Mb.min<-min(Mb.min,mean(Mij[p1,p2]))
    }}
	
    Mb <- Mb*2/(npop*(npop-1)) #mean(Mij, na.rm = T)
    #Fsts <- (Fsts - Mb)/(1 - Mb)
	for (i in 1:npop) Fi[[i]]<-2*(Fi[[i]]-Fsts[[i]])/(1-Fsts[[i]])-1
	names(Fi)<-as.character(x)
	Fis<-unlist(lapply(Fi,mean,na.rm=TRUE))
	Fis<-c(Fis,mean(Fis,na.rm=TRUE))
	
    Fsts<- c((unlist(Fsts)-Mb)/(1-Mb), mean((unlist(Fsts)-Mb)/(1-Mb), na.rm = TRUE)) #res.Mb
#	res.Mb.min<-c((Fsts-Mb.min)/(1-Mb.min),mean((Fsts-Mb.min)/(1-Mb.min),na.rm=TRUE))
#	res<-data.frame(rbind(res.Mb,res.Mb.min))
    res<-rbind(Fis,Fsts) 
    colnames(res) <- c(levels(pop), "All")
	rownames(res)<-c("Fis","Fst")
    return(list(Fi=Fi,Fs=res))
}

