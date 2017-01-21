mat2vec<-function (mat,u=FALSE) {
    n <- dim(mat)[2]
    nn <- n * (n - 1)/2
    x <- vector(length = nn)
    cum <- 0
  if (u) {
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
