test_that("diploid conversion preserves calls, missingness and names", {
  x <- matrix(c(11,12,21,22,NA,11),3,
              dimnames=list(c("a","b","c"),c("L1","L2")))
  expected <- x
  expected[] <- c(0,1,1,2,NA,0)
  expect_equal(biall2dos(x),expected)
  expect_equal(biall2dos(as.data.frame(x)),expected)
  for (z in list(matrix(NA,2,3),matrix(NA_real_,2,3))) {
    expect_identical(dim(biall2dos(z)),dim(z))
    expect_true(all(is.na(biall2dos(z))))
  }
  expect_equal(biall2dos(matrix(12,1,1)),matrix(1,1,1))
})

test_that("unsupported diploid codes are rejected rather than recoded", {
  for (bad in c(0,1,2,10,13,20,23,33,-1,11.5,Inf,-Inf)) {
    expect_error(biall2dos(matrix(c(11,bad,NA,22),2)),"encoded")
  }
})

test_that("haploid conversion validates codes and returns a matrix", {
  x <- matrix(c(1,2,NA,1),2,dimnames=list(c("a","b"),c("x","y")))
  expect_equal(biall2dos(x,diploid=FALSE),x-1)
  expect_equal(biall2dos(as.data.frame(x),diploid=FALSE),x-1)
  expect_true(all(is.na(biall2dos(matrix(NA,2,2),diploid=FALSE))))
  for (bad in c(0,3,11,12,22,-1,1.5,Inf)) {
    expect_error(biall2dos(matrix(c(1,bad),2),diploid=FALSE),"encoded")
  }
})

test_that("malformed containers and flags are rejected", {
  expect_error(biall2dos(c(11,12)),"matrix or data frame")
  expect_error(biall2dos(matrix(c("11","12"),2)),"numeric")
  for (flag in list(NA,1,"TRUE",c(TRUE,FALSE),logical()))
    expect_error(biall2dos(matrix(11,1,1),diploid=flag),"TRUE or FALSE")
})
