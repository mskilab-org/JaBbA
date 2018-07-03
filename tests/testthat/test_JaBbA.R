library(JaBbA)   
library(gUtils)
library(testthat)

context('JaBbA')

junctions = system.file("extdata", "junctions.vcf", package = 'JaBbA')
bedpe = system.file("extdata", "junctions.bedpe", package = 'JaBbA')
coverage = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')

test_that("read.junctions", {
	expect_equal(all(values(read.junctions(junctions))$tier==c(2, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 2)), TRUE)
  expect_equal(as.data.table(unlist(read.junctions(junctions))[, c()]), as.data.table(unlist(read.junctions(bedpe))[, c()]))
})

test_that("reciprocal.cycles", {
  pc = JaBbA:::reciprocal.cycles(read.junctions(junctions), paths = TRUE)
  expect_equal(unlist(pc$paths),
               structure(c(32, 33, 69, 72, 37, 38, 41, 42), names = c('631', '632', '1251', '1252', '711', '712', '751', '752')))
})

set.seed(42);
TILIM = 2
jab = JaBbA(junctions = junctions, coverage = coverage, slack = 1e4, hets = hets, tilim = TILIM, verbose = 1, overwrite = TRUE, ploidy=3.72, purity=NA)
jab.retiterate = JaBbA(junctions = junctions, coverage = coverage, slack = 1e4, tilim = TILIM, verbose = 1, overwrite = TRUE, reiterate=3, ploidy=3.72, purity=0.99)  ## reiterate > 1


test_that("JaBbA", {   
  expect_equal(length(jab$segstats), 68)
  expect_equal(round(jab$ploidy, 1), 3.6)
  expect_equal(jab$purity, 0.99)
  expect_equal(length(jab.retiterate$segstats), 72)
})


