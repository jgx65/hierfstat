fixture_ppfis <- function() {
  data.frame(pop=rep(1:2,each=4),
    L1=c(11,11,12,22,12,22,22,22),
    L2=c(12,11,12,11,22,12,22,12),
    L3=c(11,12,22,11,12,12,22,22))
}

run_ppfis <- function(dat) {
  set.seed(1)
  boot.ppfis(dat, nboot=100)$fis.ci
}

test_that("single-population intervals do not depend on labels or their type", {
  d <- fixture_ppfis()[1:4, ]
  expected <- run_ppfis(d)
  expect_true(all(is.finite(as.matrix(expected))))
  for (label in list("A", "Z", "DumPop", 20,
                     factor("Z"), factor("Z", levels=c("unused", "Z")))) {
    x <- d
    x$pop <- rep(label, nrow(x))
    result <- run_ppfis(x)
    expect_equal(unname(as.matrix(result)), unname(as.matrix(expected)))
    expect_identical(rownames(result), as.character(label))
    expect_identical(dim(result), c(1L, 2L))
  }
})

test_that("single-population genind conversion adds no output population", {
  d <- fixture_ppfis()[1:4, ]
  d$pop <- "Z"
  g <- adegenet::df2genind(d[-1], ncode=1, ploidy=2,
                          pop=factor(d$pop))
  expect_equal(run_ppfis(g), run_ppfis(d))
})

test_that("row order and unused population levels do not alter intervals", {
  d <- fixture_ppfis()
  expect_equal(run_ppfis(d[c(8,2,5,1,4,7,3,6), ]), run_ppfis(d))
  f <- d
  f$pop <- factor(f$pop, levels=c(0,1,2,3))
  expect_equal(run_ppfis(f), run_ppfis(d))
})

test_that("missing calls and a monomorphic locus can be handled", {
  d <- fixture_ppfis()[1:4, ]
  d$L1[1] <- NA
  d$mono <- 11
  a <- run_ppfis(d)
  d$pop <- "Z"
  b <- run_ppfis(d)
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)))
  expect_true(all(is.finite(as.matrix(b))))
})

test_that("a one-locus single-population result stays a data frame", {
  d <- fixture_ppfis()[1:4, 1:2]
  d$pop <- "Z"
  result <- run_ppfis(d)
  expect_identical(dim(result), c(1L, 2L))
  bs <- basic.stats(transform(d, pop=1))
  expected <- round(1-bs$Ho[1,1]/bs$Hs[1,1],4)
  expect_equal(as.numeric(result[1,]), rep(expected,2))
})
