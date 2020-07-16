"prepdata" <-
  function (data) 
  {

    if (!is.data.frame(data)) data<-as.data.frame(data) 
    nbl <- dim(data)
    nbi <-nbl[1]
    nbl<-nbl[2]
    names.data <- names(data)
    x <- matrix(numeric(nbi*nbl), ncol = nbl)
    for (i in nbl:2) {
      dumtext <- parse(text = paste("table(", paste("data[,", 
                                                    i:2, "],", sep = "", collapse = ""), "data[,1])", 
                                    sep = "", collapse = ""))
      dum <- unlist(as.vector(eval(dumtext)))
      dum1 <- dum[dum > 0]
      x[, i] <- rep(1:length(dum1), dum1)
    }
    dum <- unlist(as.vector(table(data[, 1])))
    dum1 <- dum[dum > 0]
    x[, 1] <- rep(1:length(dum1), dum1)
    x <- data.frame(x)
    names(x) <- names.data
    return(x)
  }

#################
prepdatan<-function(data){
#slower than prepdata  
  if(!is.data.frame(data)) data<-as.data.frame(data)
  nbl<-dim(data)
  nbi<-nbl[1]
  nbl<-nbl[2]
  names.data<-names(data)
  x<-matrix(numeric(nbl*nbi),ncol=nbl)
  for (i in nbl:2){
    dumtext<-parse(text=paste(paste("factor(data[,",1:(i-1),"]):",sep="",collapse=""),
                              "factor(data[,",i,"])",sep="",collapse=""))
    x[,i]<-as.integer(factor(paste(eval(dumtext))))
  }
  x[,1]<-as.integer(factor(data[,1]))
  x<-data.frame(x)
  names(x)<-names(data)
  return(x)
}

############################################
"varcompn" <-
function (data, diploid = TRUE,sorted=FALSE) 
{
    vcomp <- function(y1) {
        ss <- vector(length = nblevels)
        for (i in 1:nblevels) ss[i] <- sum(tapply(y1, ndata[, 
            i], sum)^2/table(ndata[, i]))
        temp1 <- c(sum(y1)^2/length(y1), ss)
        vec.c.ss <- temp1[2:length(temp1)] - temp1[1:(length(temp1) - 
            1)]
        meansq <- vec.c.ss/dfreed
        vc <- solve(k, meansq)
        return(vc)
    }
    
    nbf <- dim(data)[2] - 1 #nb levels
    if (!sorted){
    x <- NULL
    if (nbf > 1) 
        for (i in 1:(nbf - 1)) x <- paste(x, paste("data[,", 
            i, "],", sep = "", collapse = ""))
    no <- eval(parse(text = paste("order(", x, "data[,", nbf, 
        "])")))
    data <- data[no, ]  #for sorting the data set
    }
    y <- data[, nbf+1] # last column
    dum <- !is.na(y) # remove nas
    expl <- prepdata(cbind(data[dum, -(nbf+1)], 1:dim(data[dum, 
        ])[1])) #renumber factors to be truly nested
    if (diploid) {
        expl <- rbind(expl, expl)
        ny <- genot2al(y[dum])
        al <- 1:length(ny)
        ndata <- data.frame(expl, al, ny)
    }
    else {
        ndata <- data.frame(expl, y[dum])
    }
    rm(y)
    y <- ndata[, dim(ndata)[2]]
    nblevels <- dim(ndata)[2]-1
    names.al <- names(table(y))
    y <- as.numeric(y)
    id.al <- as.numeric(names(table(y)))
    nal <- length(id.al)
    #could use matrix rather than cbind for constructing resp
    resp<-matrix(numeric(length(y)*nal),ncol=nal)
#    resp <- as.numeric(y == id.al[1])
    for (i in 1:nal) resp[,i] <- as.numeric(y == id.al[i])
    n <- vector(length = (nblevels))
    for (i in 1:nblevels) n[i] <- max(ndata[, i])
    n <- c(1, n)
    dfreed <- n[2:(nblevels + 1)] - n[1:nblevels]
    k <- matrix(rep(0, (nblevels)^2), ncol = (nblevels))
    #could use matrix rather than cbind for constructing x
    x<-cbind(rep(1,dim(ndata)[1]),ndata[,1:nblevels])
#    x <- rep(1, length(ndata[, 1]))
#    for (i in 1:nblevels) x <- cbind(x, ndata[, i]) #OPT
    dum <- list()
    temp <- rep(1, length(y))
    for (i in 1:nblevels) dum[[i]] <- tapply(temp, x[, i], sum)
    dum[[(nblevels + 1)]] <- temp
    for (i in 2:nblevels) {
        for (j in i:nblevels) {
            temp <- length(table(x[, (i - 1)]))
            thisdum <- vector(length = 0)
            for (jj in 1:temp) thisdum <- c(thisdum, as.vector(rep(dum[[i - 
                1]][jj], length(table(x[, j][x[, (i - 1)] == 
                jj])))))
            a <- sum(dum[[j]]^2/thisdum)
            temp <- length(table(x[, i]))
            thisdum <- vector(length = 0)
            for (jj in 1:temp) thisdum <- c(thisdum, as.vector(rep(dum[[i]][jj], 
                length(table(x[, j][x[, i] == jj])))))
            b <- sum(dum[[j]]^2/thisdum)
            k[(i - 1), (j - 1)] <- (b - a)/dfreed[(i - 1)]
        }
    }
    k[, (nblevels)] <- 1
    res <- apply(resp, 2, vcomp)
    res <- data.frame(res)
    names(res) <- names.al
    res <- t(res)
    tot <- apply(res, 2, sum)
    f <- matrix(rep(0, (nblevels - 1)^2), ncol = (nblevels - 
        1))
    for (i in 1:(nblevels - 1)) {
        for (j in i:(nblevels - 1)) {
            f[i, j] <- sum(tot[i:j])/sum(tot[i:nblevels])
        }
    }
    return(list(df = dfreed, k = k, res = res, overall = tot, 
        F = f))
}

############################################
"varcompn.glob" <-
function (levels = levels, loci = loci, diploid = TRUE) 
{
    lnames <- names(loci)
    if (is.null(dim(levels))) { #only one level
      levels<-as.data.frame(levels)
        fnames <- names(levels) # "Pop"
    }
    else fnames <- names(levels)
    if (diploid) {
        fnames <- c(fnames, "Ind")  
    }
### sorting levels
    x <- NULL
    nbf<-dim(levels)[2]
    if (nbf > 1) 
      for (i in 1:(nbf - 1)) x <- paste(x, paste("levels[,", 
                                                 i, "],", sep = "", collapse = ""))
    no <- eval(parse(text = paste("order(", x, "levels[,", nbf, 
                                  "])")))
    levels <- levels[no, ]  #for sorting the data set 
    loci<-loci[no,]
    
    nloc <- dim(loci)[2] 
    nbvc<-nbf+1
    if (diploid) nbvc<-nbvc+1
    res<-matrix(numeric(nbvc*nloc),ncol=nbvc)

    for (i in 1:nloc) res[i,] <- varcompn(cbind(levels, 
        loci[, i]),diploid,sorted=TRUE)$overall 
    tot <- apply(res, 2, sum, na.rm = TRUE)
    nblevels <- length(tot)
    f <- matrix(rep(0, (nblevels - 1)^2), ncol = (nblevels - 
        1))
    for (i in 1:(nblevels - 1)) {
        for (j in i:(nblevels - 1)) {
            f[i, j] <- sum(tot[i:j])/sum(tot[i:nblevels])
        }
    }
#    fnames
    row.names(res) <- lnames
    names(tot) <- c(fnames, "Error")
    tf <- t(f)
    row.names(tf) <- fnames
    f <- t(tf)
    row.names(f) <- c("Total", fnames[-length(fnames)])
    class()
    result<-list(loc = res, overall = tot, F = f)
    class(result)<-"varcompn.glob"
    result
}

print.varcompn.glob<-function(x,...){
  print(x$F)
  invisible(x)
}
  