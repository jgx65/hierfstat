#' Reads a VCF file into a BED object
#' 
#' Reads a \url{https://samtools.github.io/hts-specs/}{Variant Call Format (VCF)} file into a BED object,
#' retaining bi-allelic SNPs only
#'
#' @usage read.VCF(fname,BiAllelic=TRUE,...)
#' 
#' @param fname VCF file name. The VCF file can be compressed (VCF.gz)
#' @param BiAllelic  Logical. If TRUE, only bi-allelic SNPs are retained,
#'    otherwise, all variant are kept
#' @param ... other arguments to pass to the function    
#'     
#'    
#' @return  A \code{\link[gaston]{bed.matrix}}   
#' 
#' @seealso \code{\link[gaston]{read.vcf}} 
#' 
#' @examples 
#'  
#'  filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
#'  x1 <- read.VCF( filepath )
#'  x1

#'
#' @export
read.VCF<-function(fname,BiAllelic=TRUE,...){
  bed<-gaston::read.vcf(fname,...)
  if(BiAllelic){
  snps<-with(bed@snps,which(nchar(A1)==1 & nchar(A2)==1))
  bed<-bed[,snps]
  }
  bed
}

