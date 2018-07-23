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

cov.gr = dt2gr(fread(coverage))
hets.gr = dt2gr(fread(hets)[, ":="(mean_high = pmax(alt, ref), mean_low = pmin(alt, ref))])
segs.gr = readRDS(segs) %$% cov.gr[, 'ratio'] %$% hets.gr[, c('mean_high', 'mean_low')]
segs.gr$mean = segs.gr$ratio
segs.gr$sd_high = segs.gr$sd_low = segs.gr$sd = 1
pp = JaBbA:::ppgrid(segs.gr, allelic = TRUE)


print(pp[1,])

test_that("ppgrid", {
  expect_equal(pp$purity[1], 0.98)
  expect_equal(pp$ploidy[1], 3.88)
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
EPGAP = 0.95
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
            epgap = EPGAP,
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
                      epgap = EPGAP)  ## reiterate > 1

test_that("karyograph", {
    kag = JaBbA:::karyograph(junctions = junc)
    expect_equal(length(kag$tile), 336)
    seqlevels(nsegs) = as.character(1:22)
    kag = JaBbA:::karyograph(junctions = junc, tile = nsegs, label.edges = TRUE)
    expect_equal(length(kag$tile), 1144)
    kag = JaBbA:::karyograph(junctions = NULL, tile = nsegs)
    expect_equal(length(kag$tile), 812)
})


list.expr = function(x)
  {
      if (is.character(x))
              paste("c('", paste(x, sep = "", collapse = "', '"), "')", sep = "")
      else
          paste("c(", paste(x, sep = "", collapse = ", "), ")", sep = "")
  }

print('jab cn')
print(list.expr(jab$segstats$cn))

print('jab junctions cn')
print(list.expr(values(jab$junctions)$cn))

print('jab purity ploidy')
print(c(jab$purity, jab$ploidy))

print('jab.reiterate cn')
print(list.expr(jab.reiterate$segstats$cn))

print('jab.reiterate junctions cn')
print(list.expr(values(jab.reiterate$junctions)$cn))

print('jab.reiterate purity ploidy')
print(c(jab.reiterate$purity, jab.reiterate$ploidy))


test_that("JaBbA", {
  expect_true( ## accounting for difference in run speed local vs travis
    identical(jab$segstats$cn,c(4, 2, 4, 3, 3, 1, 3, 28, 32, 28, 33, 28, 32, 28, 4, 24, 4, 4, 2, 4, 3, 3, 1, 3, 28, 32, 28, 33, 28, 32, 28, 4, 24, 4, 2, 1, 24, 2, 25, 2, 25, 2, 1, 24)) |
    identical(jab$segstats$cn,c(4, 5, 4, 3, 3, 1, 3, 28, 3, 4, 32, 4, 27, 4, 28, 4, 33, 4, 32, 4, 24, 4, 4, 5, 4, 3, 3, 1, 3, 28, 3, 4, 32, 4, 27, 4, 28, 4, 33, 4, 32, 4, 24, 4, 1, 1, 1, 1, 1, 1, 1, 1))
    , TRUE)

  expect_true(identical(values(jab$junctions)$cn,  c(2, 4, 5, 4, 20)),
              identical(values(jab$junctions)$cn,  c(2, 25, 28, 23, 24, 29, 28, 20)))

  expect_true(abs(jab$ploidy - 3.75)<0.01 | abs(jab$ploidy - 3.60)<0.01)

  expect_true(abs(jab$purity - .999)<0.01 | abs(jab$purity -  1.000000)<0.01)

  expect_true(identical(jab.reiterate$segstats$cn, c(3, 4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 27, 31, 40, 23, 32, 31, 32, 27, 21, 4, 3, 4, 9, 3, 3, 4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 27, 31, 40, 23, 32, 31, 32, 27, 21, 4, 3, 4, 9, 3, 1, 1, 1, 1)) |
              identical(jab.reiterate$segstats$cn, c(4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 27, 31, 40, 23, 32, 31, 32, 27, 21, 4, 3, 4, 9, 3, 4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 27, 31, 40, 23, 32, 31, 32, 27, 21, 4, 3, 4, 9, 3, 1, 1)))

  expect_true(identical(values(jab.reiterate$junctions)$cn,  c(1, 2, 1, 2, 10, 11, 4, 6, 6, 0, 9, 17, 1, 5, 1)) |
              identical(values(jab.reiterate$junctions)$cn,  c(1, 2, 1, 2, 10, 11, 4, 6, 6, 0, 9, 17, 1, 5, 1)))

  expect_true(abs(jab.reiterate$ploidy - 3.50)<0.01 |
              abs(jab.reiterate$ploidy - 3.51)<0.01)

  expect_true(abs(jab.reiterate$purity - .99)<0.01 |
              abs(jab.reiterate$purity - .99)<0.01)
  })


