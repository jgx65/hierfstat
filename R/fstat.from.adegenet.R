#################
# fstat function
#################
#
#' Wrapper for fst estimator from hierfstat package (from adegenet)
#'
#'
#' The function \code{fstat} is a wrapper for \code{varcomp.glob} of the
#' package \code{hierfstat}. For Fst, Fis and Fit, an alternative is offered by
#' \code{Fst} from the \code{pagas} package (see example).
#'
#' Let \eqn{A} and \eqn{B} be two populations of population sizes \eqn{n_A} and
#' \eqn{n_B}, with expected heterozygosity (averaged over loci) \eqn{Hs(A)} and
#' \eqn{Hs(B)}, respectively. We denote \eqn{Ht} the expected heterozygosity of
#' a population pooling \eqn{A} and \eqn{B}. Then, the pairwise \eqn{Fst}
#' between \eqn{A} and \eqn{B} is computed as:\cr
#'
#' \eqn{ Fst(A,B) = \frac{(Ht - (n_A Hs(A) + n_B Hs(B))/(n_A + n_B) )}{Ht}} \cr
#'
#' @aliases fstat FST fst 
#' @param x an object of class \linkS4class{genind}.
#' @param pop a factor giving the 'population' of each individual. If NULL, pop
#' is seeked from \code{pop(x)}. Note that the term population refers in fact
#' to any grouping of individuals'.
#' @param res.type the type of result to be returned: a \code{dist} object, or
#' a symmetric matrix
#' @param fstonly a logical stating whether only the Fst should be returned.
#' @return A vector, a matrix, or a dist object containing F statistics.
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @seealso \code{\link[adegenet]{Hs}}
#' @references Nei, M. (1973) Analysis of gene diversity in subdivided
#' populations. Proc Natl Acad Sci USA, 70: 3321-3323
#' @keywords multivariate
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @import adegenet
#'
#' @rdname fstat.from.adegenet
#'
#' @aliases fstat
#'
#' @examples
#'
#' \dontrun{
#' if(require(adegenet)){
#' data(nancycats)
#'
#'
#' ## Fst, Fis, Fit
#' ## using hierfstat
#' if(require(hierfstat)){
#' fstat(nancycats)
#' }
#'
#' ## using pegas
#' if(require(pegas)){
#' data(nancycats)
#'
#' ## conversion to pegas's format
#' as.loci(nancycats)
#'
#' ## use Fst from pegas
#' fsttab <- Fst(as.loci(nancycats))
#'
#' ## average over loci
#' apply(fsttab, 2, mean)
#' }
#' }
#'
fstat <- function(x, pop=NULL, fstonly=FALSE){
    ## cat("\nSorry, hierfstat package has been disabled - this function will be restored in a future release.\n")
    ## return(invisible())
    ## misc checks
    if(!is.genind(x)) stop("x is not a valid genind object")
    ## if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")
    if(any(x@ploidy != 2L)) stop("not implemented for non-diploid genotypes")

    if(is.null(pop)) pop <- pop(x)
    if(is.null(pop)) stop("no pop factor provided")
    if(length(pop)!=nInd(x)) stop("pop has a wrong length.")

    ## computations
    dat <- .genind2hierfstat(x)[,-1]
    res <- varcomp.glob(levels=data.frame(pop), loci=dat)$F

    if(fstonly) {res <- res[1,1]}
    return(res)
}

