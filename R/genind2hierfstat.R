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
      stop("population factor must be defined")
    } else {
      pop <- adegenet::pop(dat)
    }
  }
  
  if (dat@type!="codom") stop("data type must be codominant. Exiting")  
  ploid <- unique(dat@ploidy)
  if (length(ploid) != 1)
    stop("data must contain only diploids or only haploids. Exiting")
  if (ploid > 2L)
    stop("Data must come from diploids or haploids. Exiting")
  
  nucleotides <- c("A", "C", "G", "T")
  alleles.name <- toupper(unique(unlist(dat@all.names)))
  x <- genind2df(dat, sep = "", usepop = FALSE)
  x <- if(!all(alleles.name %in% nucleotides)) {
    do.call(cbind, lapply(x, function(x.col) as.numeric(factor(x.col))))
  } else {
    do.call(cbind, lapply(x, function(x.col) {
      x.col <- toupper(x.col)
      as.numeric(factor(x.col, levels = nucleotides))
    }))
  }
  return(data.frame(pop = pop, x))
}
  
