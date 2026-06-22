####
#'@export
"varcomp" <-
  function (data, diploid = TRUE) 
  {
    vcomp <- function(y1) {
      ss <- vector(length = nblevels)
      for (i in 1:nblevels) ss[i] <- sum(tapply(y1, ndata[, 
                                                          i], sum)^2/table(ndata[, i]))
      temp1 <- c(sum(y1)^2/length(y1), ss)
      vec.c.ss <- temp1[2:length(temp1)] - temp1[1:(length(temp1) - 
                                                      1)]
      meansq <- vec.c.ss/dfreed
      solve(k, meansq)
    }
    nbf <- dim(data)[2] - 1
    x <- NULL
    if (nbf > 1) 
      for (i in 1:(nbf - 1)) x <- paste(x, paste("data[,", 
                                                 i, "],", sep = "", collapse = ""))
    no <- eval(parse(text = paste("order(", x, "data[,", nbf, 
                                  "])")))
    data <- data[no, ]
    y <- data[, dim(data)[2]]
    dum <- !is.na(y)
    expl <- prepdata(cbind(data[dum, -dim(data)[2]], 1:dim(data[dum, 
                                                                ])[1]))
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
    nblevels <- dim(ndata)[2] - 1
    names.al <- names(table(y))
    y <- as.numeric(y)
    id.al <- as.numeric(names(table(y)))
    nal <- length(id.al)
    resp <- as.numeric(y == id.al[1])
    for (i in 2:nal) resp <- cbind(resp, as.numeric(y == id.al[i])) #OPT
    n <- vector(length = (nblevels))
    for (i in 1:nblevels) n[i] <- max(ndata[, i])
    n <- c(1, n)
    dfreed <- n[2:(nblevels + 1)] - n[1:nblevels]
    k <- matrix(rep(0, (nblevels)^2), ncol = (nblevels))
    x <- rep(1, length(ndata[, 1]))
    for (i in 1:nblevels) x <- cbind(x, ndata[, i]) #OPT
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
###
#'@export

"varcomp.glob" <-
  function (levels = levels, loci = loci, diploid = TRUE, n.cores = 1) 
{
  lnames <- names(loci)
  if (is.null(dim(levels))) {
    fnames <- "Pop"
  }
  else fnames <- names(levels)
  if (diploid) {
    fnames <- c(fnames, "Ind")
  }
  
  nloc <- dim(loci)[2]
  
  if (!is.null(n.cores) && n.cores > 1) {
    
    cores.disponiveis <- parallel::detectCores()
    
    if (n.cores >= cores.disponiveis) {
      n.cores <- max(1, cores.disponiveis - 1)
      message(paste("Adjusting n.cores to", n.cores, "to prevent system crash."))
    }
    
    message(paste("Processing", nloc, "loci using", n.cores, "cores..."))
    
    if (.Platform$OS.type == "unix") {
      lista_res <- parallel::mclapply(1:nloc, function(i) {
        varcomp(cbind(levels, loci[, i]), diploid)$overall
      }, mc.cores = n.cores)
    } else {
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl))       
      parallel::clusterExport(cl, varlist = c("levels", "loci", "diploid"), envir = environment())
      parallel::clusterEvalQ(cl, library(hierfstat))
      
      lista_res <- parallel::parLapply(cl, 1:nloc, function(i) {
        varcomp(cbind(levels, loci[, i]), diploid)$overall
      })
    }
    res <- do.call(rbind, lista_res)
    
  } else {
    lista_res <- lapply(1:nloc, function(i) {
      varcomp(cbind(levels, loci[, i]), diploid)$overall
    })
    res <- do.call(rbind, lista_res)
  }
  
  tot <- apply(res, 2, sum, na.rm = TRUE)
  nblevels <- length(tot)
  f <- matrix(rep(0, (nblevels - 1)^2), ncol = (nblevels - 1))
  for (i in 1:(nblevels - 1)) {
    for (j in i:(nblevels - 1)) {
      f[i, j] <- sum(tot[i:j])/sum(tot[i:nblevels])
    }
  }
  row.names(res) <- lnames
  names(tot) <- c(fnames, "Error")
  tf <- t(f)
  row.names(tf) <- fnames
  f <- t(tf)
  row.names(f) <- c("Total", fnames[-length(fnames)])
  return(list(loc = res, overall = tot, F = f))
}
