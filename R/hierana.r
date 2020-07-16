##################
#'@export
###################
"g.stats" <-
function (data, diploid = TRUE) 
{
    y <- data[, 2]
    x <- data[, 1]
    dum <- !is.na(y)
    if (diploid) {
        y <- genot2al(y[dum])
        x <- rep(x[dum], 2)
    }
    else {
        x <- x[dum]
        y <- y[dum]
    }
    obs <- table(x, y)
    nt <- sum(obs)
    s.r <- apply(obs, 1, sum)
    s.c <- apply(obs, 2, sum)
    expe <- s.r %*% t(s.c)/nt
    x.squared <- sum((obs - expe)^2/expe, na.rm = TRUE)
    g.stats <- 2 * sum(obs * log(obs/expe), na.rm = TRUE)
    list(obs = obs, expe = expe, x.squared = x.squared, g.stats = g.stats)
}
##################
#'@export
###################
"g.stats.glob" <-
function (data, diploid = TRUE) 
{
    dum <- (1:dim(data)[1])[complete.cases(data[, -1])]
    nl <- dim(data)[2] - 1
    g.stats.l <- vector(length = nl)
    g.stats <- 0
    for (i in 1:nl) {
        if (diploid) {
            y <- genot2al(data[dum, (i + 1)])
            x <- rep(data[dum, 1], 2)
        }
        else {
            x <- data[dum, 1]
            y <- data[dum, i + 1]
        }
        obs <- table(x, y)
        nt <- sum(obs)
        s.r <- apply(obs, 1, sum)
        s.c <- apply(obs, 2, sum)
        expe <- s.r %*% t(s.c)/nt
        g.stats.l[i] <- 2 * sum(obs * log(obs/expe), na.rm = TRUE)
    }
    g.stats <- sum(g.stats.l)
    list(g.stats.l = g.stats.l, g.stats = g.stats)
}
##################
#'@export
###################
"genot2al" <-
function (y) 
{
    if (max(y, na.rm = TRUE) <= 100) {
        modulo <- 10
    }
    else {
        if (max(y, na.rm = TRUE) <= 10000) {
            modulo <- 100
            d1<-y%/%modulo
            d2<-y%%modulo
            if (min(d1,na.rm=TRUE)>9 & max(d2,na.rm=TRUE)<10) modulo<-1000
         }
        else modulo <- 1000
    }
    al1 <- y%/%modulo
    al2 <- y%%modulo
    y.al <- c(al1, al2)
    return(y.al)
}
"prepdata" <-
function (data) 
{
#remove calls to names.data, apparently does not do anything
    nbl <- dim(data)[2]
    nbi <-dim(data)[1]
#    names.data <- names(data)
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
#    names(x) <- names.data
    return(x)
}

##################
#'@export
###################
"read.fstat.data" <-
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



