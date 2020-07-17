#' A genetic dataset from a diploid organism
#' 
#' A simple diploid dataset, with allele encoded as one digit number
#' 
#'
#a' @format A data frame with 44 rows and 6 columns:
#' \describe{
#'   \item{Pop}{Population identifier, from 1 to 6}
#'   \item{loc-1}{genotype at loc-1 (only allele 4 present)}
#'   \item{loc-2}{genotype at loc-1 (alleles 3 and 4)}
#'   \item{loc-3}{genotype at loc-1 (alleles 2, 3 and 4)}
#'   \item{loc-4}{genotype at loc-1 (alleles 1, 2, 3 and 4)}
#'   \item{loc-5}{genotype at loc-1 (only allele 4)}
#'  
#'   ...
#' }
#' 
#' @examples 
#' data(diploid)
#' basic.stats(diploid)
#' 
#' @source Given in Weir, B.S. Genetic Data Analysis. Sinauer
"diploid"
############################################################################
#' A genetic dataset from a diploid organism in a continent-island model
#' 
#' A simple diploid dataset, with allele encoded as one digit number.  
#' Up to 4 alleles per locus 
#' 
#'
#a' @format A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{Pop}{Population identifier, from 1 to 3}
#'   \item{loc.1}{genotype at loc.1}
#'   \item{loc.2}{genotype at loc.2}
#'   \item{loc.3}{genotype at loc.3}
#'   \item{loc.4}{genotype at loc.4}
#'   \item{loc.5}{genotype at loc.5}
#'  
#'   ...
#' }
#' @examples 
#' data(cont.isl)
#' allele.count(cont.isl)
#' 
#' @source generated with function sim.genot()
"cont.isl"
#########################################################################
#' A genetic dataset from a diploid organism in a continent-island model
#' 
#' A simple diploid dataset, with alleles encoded as two digits numbers.  
#' Up to 99 alleles per locus 
#' 
#'
#a' @format A data frame with 150 rows and 6 columns:
#' \describe{
#'   \item{Pop}{Population identifier, from 1 to 3}
#'   \item{loc.1}{genotype at loc.1}
#'   \item{loc.2}{genotype at loc.2}
#'   \item{loc.3}{genotype at loc.3}
#'   \item{loc.4}{genotype at loc.4}
#'   \item{loc.5}{genotype at loc.5}
#'  
#'   ...
#' }
#' 
#' @examples
#' data(cont.isl99)
#' allele.count(cont.isl99)
#' 
#' @source generated with function sim.genot(nbal=99)
"cont.isl99"
##########################################################################