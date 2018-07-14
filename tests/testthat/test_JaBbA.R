library(JaBbA)
library(gUtils)
library(testthat)

context('JaBbA')

junctions = system.file("extdata", "junctions.vcf", package = 'JaBbA')
bedpe = system.file("extdata", "junctions.bedpe", package = 'JaBbA')
coverage = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')
segs = system.file("extdata", "segs.rds", package = 'JaBbA')

test_that("read.junctions", {
    expect_equal(all(values(read.junctions(junctions))$tier==c(2, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 2)), TRUE)
    expect_equal(as.data.table(unlist(read.junctions(junctions))[, c()]), as.data.table(unlist(read.junctions(bedpe))[, c()]))
    junc.tab = fread(bedpe)[, .(chr1 = V1, pos1 = V2, chr2 = V4, pos2 = V5, str1 = V9, str2 = V10)]
    expect_equal(as.data.table(unlist(read.junctions(junc.tab))[, c()]), as.data.table(unlist(read.junctions(bedpe))[, c()]))
})

test_that("reciprocal.cycles", {
    pc = JaBbA:::reciprocal.cycles(read.junctions(junctions), paths = TRUE)
    expect_equal(unlist(pc$paths),
                 structure(c(32, 33, 69, 72, 37, 38, 41, 42),
                           names = c('631', '632', '1251', '1252', '711', '712', '751', '752')))
})

junc = read.junctions(junctions)
values(junc)$nudge = 0
junc = rep(junc, 2)

test_that("ra.merge", {
    ram = JaBbA:::ra.merge(read.junctions(junctions),
                           read.junctions(bedpe),
                           read.junctions(junctions))
    expect_equal(ncol(values(ram)), 29)
    expect_equal(length(ram), 83)
    junc2 = GenomicRanges::split(GenomicRanges::shift(unlist(junc),400),
                                 rep(c(1,2), each = length(junc)))
    ram = JaBbA:::ra.merge(junc, junc2)
    expect_equal(length(ram), 168)
})

set.seed(42);
TILIM = 3600
nsegs = readRDS(segs)
nsegs$cn = 2

jab = JaBbA(junctions = junc,
            coverage = coverage,
            seg = segs, nseg = nsegs,
            strict = TRUE,
            slack.penalty = 1e4,
            hets = hets,
            tilim = TILIM,
            cfield = 'nudge',
            verbose = 2,
            overwrite = TRUE,
            ploidy=3.72,
            purity=NA,
            epgap = 0.8,
            all.in = TRUE,
            junctions.unfiltered = junctions,
            tfield = 'nothing',
            nudge.balanced = TRUE)

hets.gr = dt2gr(fread(hets))

jab.reiterate = JaBbA(junctions = junctions,
                      coverage = coverage,
                      hets = hets.gr,
                      slack.penalty = 1e4,
                      tilim = TILIM,
                      verbose = 2,
                      overwrite = TRUE,
                      reiterate=3,
                      ploidy=3.72,
                      purity=0.99,
                      loose.penalty.mode = 'boolean',
                      epgap = 0.8)  ## reiterate > 1

test_that("karyograph", {
    kag = JaBbA:::karyograph(junctions = junc)
    expect_equal(length(kag$tile), 336)
    seqlevels(nsegs) = as.character(1:22)
    kag = JaBbA:::karyograph(junctions = junc, tile = nsegs, label.edges = TRUE)
    expect_equal(length(kag$tile), 1144)
    kag = JaBbA:::karyograph(junctions = NULL, tile = nsegs)
    expect_equal(length(kag$tile), 812)
})

test_that("JaBbA", {
  expect_equal(jab$segstats$cn, c(4, 3, 3, 1, 3, 28, 3, 32, 3, 27, 3, 27, 3, 33, 3, 32, 3, 24, 3, 4, 4, 3, 3, 1, 3, 28, 3, 32, 3, 27, 3, 27, 3, 33, 3, 32, 3, 24, 3, 4, 1, 1, 1, 1))
  expect_equal(values(jab$junctions)$cn,  c(2, 25, 29, 24, 24, 30, 29, 21))
  expect_true(abs(jab$ploidy - 3.41)<0.01)
  expect_true(abs(jab$purity - .999)<0.01)

  expect_equal(jab.reiterate$segstats$cn, c(4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 34, 27, 31, 39, 24, 32, 31, 32, 26, 19, 4, 3, 4, 10, 3, 4, 4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 34, 27, 31, 39, 24, 32, 31, 32, 26, 19, 4, 3, 4, 10, 3, 4, 1, 1, 1, 1))
  expect_equal(values(jab.reiterate$junctions)$cn,  c(1, 2, 1, 2, 10, 11, 4, 7, 7, 0, 8, 15, 1, 6, 1))
  expect_true(abs(jab.reiterate$ploidy - 3.70)<0.01)
  expect_true(abs(jab.reiterate$purity - .99)<0.01)
})


