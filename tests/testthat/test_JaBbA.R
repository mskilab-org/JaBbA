library(JaBbA)
library(gUtils)
library(testthat)

context('JaBbA')

junctions = system.file("extdata", "junctions.vcf", package = 'JaBbA')
bedpe = system.file("extdata", "junctions.bedpe", package = 'JaBbA')
coverage = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')
segs = system.file("extdata", "segs.rds", package = 'JaBbA')

system('rm -rf JaBbA')

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
TILIM = 3600
jab = JaBbA(junctions = junctions, coverage = coverage, seg = segs, slack.penalty = 1e4, hets = hets, tilim = TILIM, verbose = 2, overwrite = TRUE, ploidy=3.72, purity=NA, epgap = 0.1)
#jab2 = readRDS('JaBbA/jabba.rds')
jab.reiterate = JaBbA(junctions = junctions, coverage = coverage, seg = segs, slack.penalty = 1e4, tilim = TILIM, verbose = 2, overwrite = TRUE, reiterate=3, ploidy=3.72, purity=0.99, epgap = 0.1)  ## reiterate > 1
#jab.reiterate2 = readRDS('JaBbA/jabba.rds')

print(length(jab$segstats))
print(jab$ploidy)
print(jab$purity)
print(length(jab.reiterate$segstats))


test_that("JaBbA", {
  expect_equal(length(jab$segstats), 68)
  expect_equal(round(jab$ploidy, 1), 3.6)
  expect_equal(jab$purity, 0.99)
  expect_equal(length(jab.reiterate$segstats), 68)
})


