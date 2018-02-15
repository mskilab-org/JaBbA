library(JaBbA)

junctions = system.file("extdata", "junctions.vcf", package = 'JaBbA')
coverage = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')

test_that("JaBbA", {
  set.seed(42);
  jab = JaBbA(junctions = junctions, coverage = coverage, tilim = 10, verbose = 1, overwrite = TRUE)
  expect_equal(all(jab$segstats$cn == c(2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1, 2, 8, 15, 8, 13, 18, 19, 14, 19, 25, 13, 19, 15, 20, 16, 15, 3, 2, 3, 7, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1, 2, 8, 15, 8, 13, 18, 19, 14, 19, 25, 13, 19, 15, 20, 16, 15, 3, 2, 3, 7, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)), TRUE)
  expect_equal(abs(jab$ploidy-2.550274)<1e-6, TRUE)
  expect_equal(abs(jab$purity-1)<1e-6, TRUE)
    
  set.seed(42); jab = JaBbA(junctions = junctions, coverage = coverage, hets = hets, tilim = 10, verbose = 1, overwrite = TRUE)
  expect_equal(all(jab$segstats$cn == c(3, 4, 3, 4, 3, 4, 5, 4, 5, 2, 5, 4, 3, 2, 3, 4, 3, 1, 3, 11, 22, 11, 18, 25, 27, 20, 27, 35, 19, 27, 21, 29, 23, 21, 5, 4, 3, 4, 10, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 5, 4, 5, 2, 5, 4, 3, 2, 3, 4, 3, 1, 3, 11, 22, 11, 18, 25, 27, 20, 27, 35, 19, 27, 21, 29, 23, 21, 5, 4, 3, 4, 10, 3, 4, 3, 4, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)), TRUE)
  expect_equal(abs(jab$ploidy-3.732664)<1e-6, TRUE)
  expect_equal(abs(jab$purity-0.9949749)<1e-6, TRUE)
})

test_that("read.junctions", {
  expect_equal(all(values(read.junctions(junctions))$tier==c(2, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 2)), TRUE)
})



