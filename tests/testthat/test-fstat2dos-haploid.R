test_that("haploid alleles produce one copy rather than two", {
  x <- matrix(c(1,2,1,2,2,1),3,
              dimnames=list(c("a","b","c"),c("L1","L2")))
  expected <- cbind(c(1,0,1),c(0,1,0),c(0,0,1),c(1,1,0))
  dimnames(expected) <- list(rownames(x),c("L1.1","L1.2","L2.1","L2.2"))
  expect_equal(fstat2dos(x,diploid=FALSE),expected)
  expect_equal(fstat2dos(as.data.frame(x),diploid=FALSE),expected)
  expect_equal(rowSums(fstat2dos(x,diploid=FALSE)),setNames(rep(2,3),rownames(x)))
})

test_that("haploid missing calls remain missing for every allele", {
  x <- matrix(c(3,NA,1,3),4,dimnames=list(NULL,"L"))
  y <- fstat2dos(x,diploid=FALSE)
  expect_identical(colnames(y),c("L.1","L.3"))
  expect_true(all(is.na(y[2,])))
  expect_equal(rowSums(y[-2,,drop=FALSE]),rep(1,3))
  expect_equal(dim(y),c(4L,2L))
})

test_that("monomorphic and entirely missing haploid loci retain matrix shape", {
  x <- data.frame(mono=c(7,7,NA),empty=rep(NA,3))
  y <- fstat2dos(x,diploid=FALSE)
  expect_identical(dim(y),c(3L,1L))
  expect_identical(colnames(y),"mono.7")
  expect_equal(as.numeric(y),c(1,1,NA))
  expect_identical(dim(fstat2dos(matrix(NA,3,2),diploid=FALSE)),c(3L,0L))
  expect_identical(dim(fstat2dos(matrix(1,1,1),diploid=FALSE)),c(1L,1L))
  expect_identical(colnames(fstat2dos(matrix(c(9,2,9),3),diploid=FALSE)),c("l.1.2","l.1.9"))
})

test_that("unsupported haploid identifiers are rejected explicitly", {
  for (bad in c(0,-1,1.5,Inf,-Inf))
    expect_error(fstat2dos(matrix(c(1,bad),2),diploid=FALSE),"positive integers")
  expect_error(fstat2dos(matrix(c("A","B"),2),diploid=FALSE),"positive integers")
})
