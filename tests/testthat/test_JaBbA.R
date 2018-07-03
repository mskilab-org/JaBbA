library(JaBbA)   
library(gUtils)
library(testthat)
## library(Rcplex)
context("test JaBbA functionality\n")

junctions = system.file("extdata", "junctions.vcf", package = 'JaBbA')
coverage = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')


test_that("read.junctions", {
	expect_equal(all(values(read.junctions(junctions))$tier==c(2, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 2)), TRUE)
})


set.seed(42);

jab = JaBbA(junctions = junctions, coverage = coverage, tilim = 10, verbose = 1, overwrite = TRUE, ploidy=3.72, purity=0.99)
jab_chromoplexy = JaBbA(junctions = junctions, coverage = coverage, tilim = 10, verbose = 1, overwrite = TRUE, nudge.balanced=TRUE, ploidy=3.72, purity=0.99)
jab2 = JaBbA(junctions = junctions, coverage = coverage, hets = hets, tilim = 10, verbose = 1, overwrite = TRUE, ploidy=3.72, purity=0.99)
jab.retiterate = JaBbA(junctions = junctions, coverage = coverage, tilim = 10, verbose = 1, overwrite = TRUE, reiterate=3, ploidy=3.72, purity=0.99)  ## reiterate > 1

test_that("JaBbA", {
    
    expect_equal(length(jab$segstats), 96)
    expect_equal(round(jab$ploidy, 1), 3.7)
    expect_equal(jab$purity, 0.99)

    expect_equal(length(jab_chromoplexy$segstats), 100)
    expect_equal(round(jab_chromoplexy$ploidy, 1), 3.8)
    expect_equal(jab_chromoplexy$purity, 0.99)

    expect_equal(length(jab2$segstats), 96)
    expect_equal(round(jab2$ploidy, 1), 3.7)

    expect_equal(length(jab.retiterate$segstats), 100)

})





