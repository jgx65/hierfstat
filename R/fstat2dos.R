#####################################################
#'
#' Converts a hierfstat genetic data frame to dosage data
#' 
#' Converts a hierfstat genetic data frame to dosage. For each allele at each locus, 
#' allelic dosage (number of copies of the allele) is reported. The column name is the allele 
#' identifier 
#' 
#' 
#' @usage fstat2dos(dat,diploid=TRUE)
#' 
#' @param dat data frame with genetic data without the first column (population identifier)
#' @param diploid whether the data set is from a diploid organism
#' 
#' @return a matrix with \eqn{\sum_l n_l^a} columns (where \eqn{n_l^a} is the number of alleles 
#' at locus l), as many rows as individuals, and containing the number of copies (dosage) of the 
#' corresponding allele
#' 
#' @examples
#' \dontrun{
#' dat<-sim.genot(nbal=5,nbloc=10)
#' dos<-fstat2dos(dat[,-1])
#' dim(dos) 
#' wc(dat)
#' fst.dosage(dos,pop=dat[,1])
#' 
#' } 
#' @export
#' 
#####################################################

fstat2dos<-function (dat, diploid = TRUE) 
{
  if (diploid) 
    dat <- getal.b(dat)
  lnames <- colnames(dat)
  n.loc <- dim(dat)[2]
  #to avoid getting an array in al.list rather than a list, 
  #when all loci have the same number if alleles
  all.count.p.loc<-apply(dat, 2, table)
  if (!is.list(all.count.p.loc)){
      tmp<-dim(dat)
      tmp[2]<-tmp[2]+1
      dat.tmp<-array(data=0,dim=tmp)
      dat.tmp[,-(n.loc+1),]<-dat
      all.count.p.loc<-apply(dat.tmp,2,table)[-(n.loc+1)]
      rm(dat.tmp)
  }
  
  al.list <- lapply(all.count.p.loc, function(x) as.integer(names(x)))
  nb.al <- cumsum(unlist(lapply(al.list, length)))
  tally <- function(x) sapply(al.list[[x]], function(y) rowSums(dat[, 
    x, ] == y))
  res <- lapply(1:n.loc, tally)
  if (!is.null(lnames)) {
    lnames <- unlist(lapply(1:n.loc, function(x) paste(lnames[[x]], 
      al.list[[x]], sep = ".")))
  }
  else {
    lnames <- unlist(lapply(1:n.loc, function(x) paste("l", x, 
      al.list[[x]], sep = ".")))
  }
  allres<-matrix(unlist(res),nrow=dim(dat)[1])
  colnames(allres) <- lnames
  as.matrix(allres)
}


#################################################
#'
#'Converts bi-allelic SNPs from hierfstat format to dosage format
#' 
#' Converts bi-allelic SNPs hierfstat format to dosage format, the number of alternate allele copies at a locus
#' for an individual, i.e. 11 -> 0; 12 or 21 >1 and 22 ->2
#' 
#' @usage biall2dos(dat,diploid=TRUE)
#' 
#' @param dat a hierfstat data frame without the first column (the population identifier), 
#' individuals in rows, columns with individual genotypes encoded as 11, 12, 21 and 22.
#' A numeric matrix is also accepted. For haploid data, use allele codes 1 and 2.
#' Missing calls must be NA; zero, partial calls and other codes are rejected.
#' @param diploid whether the data set is from a diploid organism
#' 
#' @return a matrix containing allelic dosages (0, 1 or 2 for diploids;
#' 0 or 1 for haploids), preserving missing values and dimension names
#' 
#' @examples
#' \dontrun{
#' biall2dos(sim.genot(nbal=2,nbloc=10)[,-1])  # a 10 column matrix
#' }
#' 
#' @export
#' 
#########################################################

biall2dos<-function(dat,diploid=TRUE){
  if (!is.logical(diploid) || length(diploid)!=1L || is.na(diploid))
    stop("diploid must be TRUE or FALSE")
  if (!is.matrix(dat) && !is.data.frame(dat))
    stop("dat must be a numeric matrix or data frame")
  dat <- as.matrix(dat)
  if (!is.numeric(dat) && !(is.logical(dat) && all(is.na(dat))))
    stop("dat must contain numeric genotype codes")
  allowed <- if (diploid) c(11,12,21,22) else c(1,2)
  if (any(!is.na(dat) & !(dat %in% allowed)))
    stop(if (diploid) "genotypes must be encoded as 11, 12, 21 or 22"
         else "haploid alleles must be encoded as 1 or 2")
  if (!diploid) return(dat-1)
  return(dat%/%10+dat%%10-2)
}
