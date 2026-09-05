test_that("mixed widths are decoded independently by locus", {
  x<-data.frame(L1=c(150150,150142,134134,150134),
                L2=c(8882,8882,8880,8882),L3=c(11,12,22,21))
  y<-getal.b(x)
  expect_equal(y[,1,1],c(150,150,134,150))
  expect_equal(y[,1,2],c(150,142,134,134))
  expect_equal(y[,2,1],rep(88,4))
  expect_equal(y[,2,2],c(82,82,80,82))
  expect_equal(y[,3,1],c(1,1,2,2))
  expect_equal(y[,3,2],c(1,2,2,1))
  expect_equal(getal.b(x[c(3,1,2)])[,c(2,3,1),],y)
  for (j in seq_len(ncol(x)))
    expect_equal(getal.b(x[,j,drop=FALSE])[,1,],y[,j,])
})
test_that("missing loci and calls retain their array positions", {
  x<-data.frame(empty=c(NA,NA,NA),typed=c(150150,NA,134150))
  y<-getal.b(x)
  expect_identical(dim(y),c(3L,2L,2L))
  expect_true(all(is.na(y[,1,])))
  expect_true(all(is.na(y[2,2,])))
  expect_equal(y[c(1,3),2,1],c(150,134))
  expect_true(all(is.na(getal.b(matrix(NA,2,2)))))
})
test_that("single-row and ordinary packed genotypes decode", {
  expect_identical(dim(getal.b(data.frame(L=1122))),c(1L,1L,2L))
  expect_equal(as.numeric(getal.b(data.frame(L=1122))),c(11,22))
  expect_equal(as.numeric(getal.b(data.frame(L=c(1001,2002)))),c(1,2,1,2))
  expect_error(getal.b(data.frame(L=1000000)),"3 digits")
})

test_that("allele counts and wc use the same decoding as getal.b", {
  d<-data.frame(pop=c(1,1,2,2),
                L1=c(150150,150142,134134,150134),
                L2=c(8882,8882,8880,8882))
  expected<-d
  expected$L2<-c(88082,88082,88080,88082)
  expect_equal(getal(d),getal(expected))
  expect_equal(allele.count(d),allele.count(expected))
  expect_equal(wc(d)$FST,wc(expected)$FST)
  expect_equal(wc(d)$FIS,wc(expected)$FIS)
  expect_equal(basic.stats(d)$Ho,basic.stats(expected)$Ho)
})
