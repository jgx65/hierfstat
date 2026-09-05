test_that("write.ped writes physical positions without changing genotypes", {
  d <- data.frame(pop=c(1,1,2), A=c(11,12,22), B=c(12,22,11))
  a <- tempfile()
  b <- tempfile()
  on.exit(unlink(paste0(rep(c(a,b), each=2), c(".map",".ped"))))
  write.ped(d, fname=a)
  write.ped(d, fname=b, loc.pos=c(100,200000000))
  m <- read.table(paste0(b,".map"))
  expect_equal(m[[4]], c(100,200000000))
  expect_equal(m[[2]], c("LA","LB"))
  expect_true(all(m[[1]] == 0 & m[[3]] == 0))
  expect_equal(read.table(paste0(a,".map"))[[4]], c(0,0))
  expect_identical(readLines(paste0(a,".ped")), readLines(paste0(b,".ped")))
  write.ped(d[,-1], pop=d$pop, fname=b, loc.pos=c(100,200000000))
  expect_equal(read.table(paste0(b,".map")), m)
  write.ped(d[,1:2], fname=b, loc.pos=42)
  expect_equal(read.table(paste0(b,".map"))[[4]], 42)
  expect_equal(ncol(read.table(paste0(b,".ped"))), 8L)
})

test_that("invalid positions are rejected before output is written", {
  d <- data.frame(pop=c(1,1), A=c(11,12), B=c(12,22))
  for (pos in list(1, numeric(), c(1,2,3), c(NA,2), c(Inf,2),
                   c(-1,2), c(1.5,2), c("1","2"), c(TRUE,FALSE),
                   matrix(c(1,2),1))) {
    f <- tempfile()
    expect_error(write.ped(d, fname=f, loc.pos=pos), "loc.pos")
    expect_false(file.exists(paste0(f,".map")))
    expect_false(file.exists(paste0(f,".ped")))
  }
})
