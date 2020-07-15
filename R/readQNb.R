##################
#'
#'
#'@export
###################
qn2.read.fstat<-function (fname, na.s = c("NA", "NaN")) {
  #reads quantinemo extended format
  #fname is the file name of quantinemo output
  #returns a list with
  #dat: a hierfstat data frame
  #sex: the sex of the individuals
  #ped: the individuals pedigree
  #W: the individual fitness
  
  x <- scan(fname, n = 4)
  nloc <- x[2]
  lnames <- scan(fname, what = character(), skip = 1, nlines = nloc)
  lnames <- c("Pop", lnames)
  dat <- scan(fname, skip = nloc + 1, na.strings = na.s, comment.char = "_")
  dat <- data.frame(matrix(dat, ncol = nloc + 4, byrow = TRUE))
  age <- dat[, nloc + 2]
  sex <- dat[age == 2, nloc + 3]
  asex <- character(length(sex))
  asex[sex == 0] <- "M"
  asex[sex == 1] <- "F"
  dat <- dat[age == 2, 1:(nloc + 1)]
  names(dat) <- lnames
  x<-readLines(fname)
  x<-x[-c(1:(nloc+1))]
  tmp<-matrix(unlist(lapply(x,
      function(y) unlist(strsplit(y,split=" "))[(nloc+4):(nloc+6)])),ncol=3,byrow=TRUE)
  ped<-data.frame(ind=as.character(tmp[,1]),dam=as.character(tmp[,2]),sire=as.character(tmp[,3]),stringsAsFactors = FALSE)
  ped<-ped[age==2,]
  w<-as.numeric(unlist(lapply(x,function(y) unlist(strsplit(y,split=" "))[nloc+7])))
  w<-w[age==2]
  return(list(dat = dat, sex = asex,ped=ped,W=w))
}
