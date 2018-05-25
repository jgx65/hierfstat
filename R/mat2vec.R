#############################################################################################
#' @title Creates a vector from a matrix
#'
#' @description creates a vector from a matrix 
#'  
#' @usage vec2mat(x,upper=FALSE)
#'
#' @param x a diagonal matrix
#' @param upper whether the upper triangular matrix is to be copied to the vector
#'
#' @return a vector 
#'
#' @examples{
#' 
#'  mat2vec(matrix(1:16,nrow=4))
#'  mat2vec(matrix(1:16,nrow=4),upper=TRUE)
#'}
#' @export

#################################################################
mat2vec<-function (mat,upper=FALSE) {
    n <- dim(mat)[2]
    nn <- n * (n - 1)/2
    x <- vector(length = nn)
    cum <- 0
  if (upper) {
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        cum <- cum + 1
        x[cum] <- mat[i, j]
      }
    }
  }
  else {
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      cum <- cum + 1
      x[cum] <- mat[i, j]
    }
  }
  }
  return(x)
}
