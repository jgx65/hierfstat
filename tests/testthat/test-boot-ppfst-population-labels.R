ppfst_label_fixture <- function() {
  data.frame(pop=rep(1:3,each=4),
    L1=c(11,11,12,22,12,22,22,22,11,12,11,12),
    L2=c(12,11,12,11,22,12,22,12,11,11,12,12),
    L3=c(11,12,22,11,12,12,22,22,11,12,12,22))
}
ppfst_label_run <- function(d) {
  set.seed(87)
  boot.ppfst(d,nboot=30)
}

test_that("population labels do not change numerical intervals", {
  d <- ppfst_label_fixture()
  expected <- ppfst_label_run(d)
  expect_true(all(is.finite(expected$ll[upper.tri(expected$ll)])))
  for (labels in list(c(10,20,30), c("A","B","C"),
                     factor(c("A","B","C")),
                     factor(c("A","B","C"),levels=c("unused","A","B","C","empty")))) {
    x <- d
    x$pop <- labels[x$pop]
    actual <- ppfst_label_run(x)
    expect_equal(unname(actual$ll),unname(expected$ll))
    expect_equal(unname(actual$ul),unname(expected$ul))
    expect_identical(rownames(actual$ll),as.character(labels))
    expect_identical(colnames(actual$ll),as.character(labels))
  }
})

test_that("input row order does not change the population mapping", {
  d <- ppfst_label_fixture()
  d$pop <- c("A","B","C")[d$pop]
  a <- ppfst_label_run(d)
  b <- ppfst_label_run(d[c(12,2,7,1,8,4,10,5,3,11,6,9),])
  expect_equal(a$ll,b$ll)
  expect_equal(a$ul,b$ul)
})

test_that("factor level ordering is retained without empty populations", {
  d <- ppfst_label_fixture()
  d$pop <- factor(c("A","B","C")[d$pop],levels=c("C","unused","A","B"))
  a <- ppfst_label_run(d)
  numeric <- d
  numeric$pop <- match(as.character(d$pop),c("C","A","B"))
  b <- ppfst_label_run(numeric)
  expect_identical(rownames(a$ll),c("C","A","B"))
  expect_equal(unname(a$ll),unname(b$ll))
  expect_equal(unname(a$ul),unname(b$ul))
})

test_that("missing identifiers and single observed populations fail clearly", {
  d <- ppfst_label_fixture()
  d$pop[1] <- NA
  expect_error(ppfst_label_run(d),"non-missing")
  d <- ppfst_label_fixture()[1:4,]
  d$pop <- factor(d$pop,levels=1:3)
  expect_error(ppfst_label_run(d),"at least two")
})
