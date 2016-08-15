########################################################
#' @title Converts genind objects from adegenet into a hierfstat data frame
#' @description  Converts genind objects from adegenet into a hierfstat data frame
#' @usage genind2hierfstat(dat,pop=NULL)
#' @param dat a genind object
#' @param pop a vector containing the population to which each individual belongs. 
#' If pop=NULL, pop taken from slot pop of the genind object
#'  
#' @return a data frame with nloci+1 columns and ninds rows.  The first column
#' contains the population identifier, the following the genotypes at each locus
#' 
#' @examples
#' 
#' \dontrun{
#' library(adegenet)
#' data(nancycats)
#' genind2hierfstat(nancycats)
#' basic.stats(nancycats)
#' genet.dist(nancycats)
#' data(H3N2)
#' basic.stats(genind2hierfstat(H3N2,pop=rep(1,dim(H3N2@@tab)[1])),diploid=FALSE)
#' }
#' 
#' @export
##########################################################
genind2hierfstat<-function(dat,pop=NULL){
  if (!is.genind(dat)) stop("dat must be a genind object. Exiting")
  
  if(is.null(pop)){
    if (is.null(adegenet::pop(dat))){
      stop("population factor must be defined") #warning instead of stop?
    } else {
      pop <- adegenet::pop(dat)
    }
  }
  
  if (dat@type!="codom") stop("data type must be codominant. Exiting")  
  ploid<-unique(dat@ploidy)
  if (length(ploid)!=1) stop("data must contain only diploids or only haploids. Exiting")
  if (ploid>2L) stop("Data must come from diploids or haploids. Exiting")
  alleles.name<-toupper(unique(unlist(adegenet::alleles(dat))))
  ids <- adegenet::indNames(dat)
  
  # convert alleles
  x <- if(!all(alleles.name %in% c("A", "C", "G", "T"))) {
    if(length(grep("[[:alpha:]|[:punct:]]", alleles.name)) > 0) {
      # for loci where alleles are not all numeric nor all nucleotides
      max.length <- max(sapply(alleles(dat), length))
      digits <- floor(log10(max.length))
      dat <- as.matrix(genind2df(dat, sep = "", usepop = FALSE,
                                 oneColPerAll = TRUE))
      dat <- apply(dat, 2, function(a) ifelse(a == "NA", NA, a))
      # cycle through each locus
      do.call(cbind, lapply(seq(ploid, ncol(dat), by = ploid),
            function(end) {
             start <- end - ploid + 1
             allelesid <- sort(unique(as.vector(dat[, start:end])))
             # match alleles of each genotype with sorted order: index is new  allele
             apply(dat[, start:end, drop = FALSE], 1, function(gntp) {
                  if(any(is.na(gntp))) return(NA)
                  gntp <- match(gntp, allelesid)
                  # zero pad numeric alleles and return collapsed genotype
                  gntp <- formatC(gntp, digits = digits, flag = "0", mode ="integer")
                  as.integer(paste(gntp, collapse = ""))
                  })
            }))
    }
    else { # all alleles are numeric
      # check if all alleles code have same number of digits
      tmp1<-unique(nchar(unlist(dat@all.names)))
      if (length(tmp1)>1){
        dig<-max(tmp1)
        allnames<-lapply(dat@all.names,formatC,width=dig,flag="0",mode="integer")
        dat@all.names<-allnames
      }  
      dat <- genind2df(dat, sep = "", usepop = FALSE)
      do.call(cbind, lapply(dat, as.integer))
    }
  } else { # all alleles are nucleotides
    dat <- genind2df(dat, sep = "", usepop = FALSE)
    do.call(cbind, lapply(dat, function(a) {
      a <- gsub("[aA]", "1", a)
      a <- gsub("[cC]", "2", a)
      a <- gsub("[gG]", "3", a)
      a <- gsub("[tT]", "4", a)
      as.integer(a)
    }))
  }
  x<-data.frame(pop=pop,x)
  rownames(x) <- ids
  return(x)
}
