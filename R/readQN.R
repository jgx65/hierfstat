#################################
#'Read QuantiNemo extended format for genotype files
#'  
#'Read QuantiNemo (\url{http://www2.unil.ch/popgen/softwares/quantinemo/}) genotype files extended format (option 2)
#'  
#' @usage qn2.read.fstat(fname, na.s = c("NA","NaN"))
#' @param fname quantinemo file name
#' @param na.s  na string used
#' @return dat  a data frame with nloc+1 columns, the first being the population
#'  to which the individual belongs and the next being the genotypes, one column per locus; 
#'  and ninds rows
#' @return sex  the sex of the individuals
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#' @seealso \code{\link{read.fstat}}
#' @references 
#' \href{https://pubmed.ncbi.nlm.nih.gov/30816926/}{Neuenschwander S, Michaud F, Goudet J (2019)} 
#' QuantiNemo 2: a Swiss knife to simulate complex demographic and genetic scenarios, 
#' forward and backward in time. Bioinformatics 35:886
#' 
#' \href{https://academic.oup.com/bioinformatics/article/24/13/1552/237901}{Neuenschwander S, Hospital F, Guillaume F, Goudet J (2008)} 
#' quantiNEMO: an individual-based program to simulate quantitative traits with explicit 
#' genetic architecture in a dynamic metapopulation. Bioinformatics 24:1552
#'
#' @examples 
#'   dat<-qn2.read.fstat(system.file("extdata","qn2_sex.dat",package="hierfstat"))
#'   sexbias.test(dat[[1]],sex=dat[[2]])
#' @export
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
