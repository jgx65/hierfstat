#' @title Fills a triangular matrix from the inputed vector
#'
#' @description Fills a triangular matrix from the inputed vector 
#'  
#' @usage vec2mat(x,diag=FALSE,upper=FALSE)
#'
#' @param x a vector
#' @param diag whether the vector contains the diagonal elements
#' @param upper whether the vector contains the upper trinagular matrix elements
#'
#' @return a matrix 
#'
#' @examples{
#'  vec2mat(1:10)
#'  vec2mat(1:10,diag=TRUE)
#'  vec2mat(1:10,upper=TRUE)
#'}
#' @export

#################################################################
vec2mat<-function(x,diag=FALSE,upper=FALSE){
	#fill a lower triangular matrix from a vector and copy it to upper triangle
nn<-length(x)
if (diagT) n<-(-1+(1+8*nn)^0.5)/2 #dim of the matrix
else n<-(1+(1+8*nn)^0.5)/2 #dim of the matrix
mat<-matrix(rep(0,n*n),ncol=n,nrow=n)
cum<-0
if (!upper){
if(diagT){
for (i in 2:n) {
	for (j in 1:i){
		cum<-cum+1
		mat[i,j]<-x[cum]
		mat[j,i]<-mat[i,j]
	}
}
} else
{
for (i in 2:n) {
	for (j in 1:(i-1)){
		cum<-cum+1
		mat[i,j]<-x[cum]
		mat[j,i]<-mat[i,j]
	}
}
}
}
else{
if(diagT){
for (i in 1:n) {
	for (j in i:n){
		cum<-cum+1
		mat[i,j]<-x[cum]
		mat[j,i]<-mat[i,j]
	}
}
} else
{
for (i in 1:(n-1)) {
	for (j in (i+1):n){
		cum<-cum+1
		mat[i,j]<-x[cum]
		mat[j,i]<-mat[i,j]
	}
}
}
}
return(mat)
}
