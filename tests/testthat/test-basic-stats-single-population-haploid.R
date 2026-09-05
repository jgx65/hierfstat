test_that("single-population haploid statistics retain diversity and counts", {
  d <- data.frame(pop=rep("A",4), L1=c(1,1,2,2), L2=c(1,2,3,3))
  r <- basic.stats(d, diploid=FALSE)
  expect_equal(as.numeric(r$Hs), c(0.6667,0.8333))
  expect_equal(as.numeric(r$n.ind.samp), c(4,4))
  expect_identical(r$Ho, NA_real_)
  expect_identical(r$Fis, NA_real_)
  expect_true(all(is.na(r$perloc$Ho)))
  expect_true(all(is.na(r$perloc$Fis)))
  expect_equal(colnames(r$Hs), "A")
  for (label in list(rep("Z",4), rep(10,4), rep("DumPop",4),
                     factor(rep("A",4),levels=c("A","unused")))) {
    d$pop <- label
    s <- basic.stats(d, diploid=FALSE)
    expect_equal(as.numeric(s$Hs), as.numeric(r$Hs))
    expect_equal(as.numeric(s$n.ind.samp), as.numeric(r$n.ind.samp))
    expect_equal(colnames(s$Hs), unique(as.character(label)))
    expect_identical(s$Fis, NA_real_)
  }
})

test_that("single-population diploid factor and dummy-name labels work", {
  d <- data.frame(pop=rep("A",4), L1=c(11,12,22,12), L2=c(11,11,12,22))
  r <- basic.stats(d)
  for (label in list(factor(rep("A",4),levels=c("A","unused")),
                     rep("DumPop",4))) {
    d$pop <- label
    s <- basic.stats(d)
    expect_equal(as.numeric(s$Hs),as.numeric(r$Hs))
    expect_equal(as.numeric(s$Fis),as.numeric(r$Fis))
    expect_equal(colnames(s$Hs),unique(as.character(label)))
  }
})
