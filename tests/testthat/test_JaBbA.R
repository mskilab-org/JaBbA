library(JaBbA)
library(gUtils)
library(gTrack)
library(testthat)

context('JaBbA')

print(sessionInfo())

## load testing data
juncs.fn = system.file("extdata", "junctions.vcf", package = 'JaBbA')
bedpe = system.file("extdata", "junctions.bedpe", package = 'JaBbA')
cov.fn = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')
segs = system.file("extdata", "segs.rds", package = 'JaBbA')
blacklist.junctions = system.file("extdata", "blacklist.junctions.rds", package = 'JaBbA')
whitelist.junctions = system.file("extdata", "whitelist.junctions.rds", package = 'JaBbA')
blacklist.coverage = system.file("extdata", "hg19.blacklist.coverage.rds", package = 'JaBbA')

Sys.setenv(DEFAULT_GENOME = paste0(system.file(
               "extdata",
               package = 'JaBbA'),
               "/human_g1k_v37.chrom.sizes")) ## must do this           
print(Sys.getenv("DEFAULT_GENOME"))
## print("This is the default seqlengths: ")
## print(hg_seqlengths())
## print("IS IT A VECTOR NOW??????!!!!!!")
## print(is.vector(hg_seqlengths()))

test_that("read.junctions", {
    expect_equal(all(values(read.junctions(juncs.fn))$tier==c(2, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 2)), TRUE)
    expect_equal(as.data.table(unlist(read.junctions(juncs.fn))[, c()]), as.data.table(unlist(read.junctions(bedpe))[, c()]))
    junc.tab = fread(bedpe)[, .(chr1 = V1, pos1 = V2, chr2 = V4, pos2 = V5, str1 = V9, str2 = V10)]
    expect_equal(as.data.table(unlist(read.junctions(junc.tab))[, c()]), as.data.table(unlist(read.junctions(bedpe))[, c()]))
    
})

test_that("reciprocal.cycles", {
    pc = JaBbA:::reciprocal.cycles(read.junctions(juncs.fn), paths = TRUE)
    expect_equal(unlist(pc$paths),
                 structure(c(32, 33, 69, 72, 37, 38, 41, 42),
                           names = c('631', '632', '1251', '1252', '711', '712', '751', '752')))
})

cov.gr = dt2gr(fread(cov.fn))
cov.gr$read_depth = cov.gr$ratio
hets.gr = dt2gr(fread(hets)[, ":="(mean_high = pmax(alt.count.t, ref.count.t), mean_low = pmin(alt.count.t, ref.count.t))])
segs.gr = readRDS(segs) %$% cov.gr[, 'ratio'] %$% hets.gr[, c('mean_high', 'mean_low')]
segs.gr$mean = segs.gr$ratio
segs.gr$sd_high = segs.gr$sd_low = segs.gr$sd = 1

pp = JaBbA:::ppgrid(segs.gr, allelic = TRUE)
## print(pp[1,])

test_that("ppgrid", {
    expect_equal(pp$purity[1], 0.98)
    expect_equal(pp$ploidy[1], 3.88)
})

juncs = read.junctions(juncs.fn)
values(juncs)$nudge = 0
junc = rep(juncs, 2)

test_that("ra.merge", {
    ram = JaBbA:::ra.merge(read.junctions(juncs.fn),
                           read.junctions(bedpe),
                           read.junctions(juncs.fn))
    expect_equal(ncol(values(ram)), 29)
    expect_equal(length(ram), 83)
    junc2 = GenomicRanges::split(GenomicRanges::shift(unlist(junc),400),
                                 rep(c(1,2), each = length(junc)))
    ram = JaBbA:::ra.merge(junc, junc2)
    expect_equal(length(ram), 168)
})

set.seed(42);
TILIM = 60
EPGAP = 0.95

nsegs = readRDS(segs)
nsegs$cn = 2

hets.gr = dt2gr(fread(hets))

test_that("karyograph", {
    kag = JaBbA:::karyograph(junctions = junc)
    expect_equal(length(kag$tile), 336)
    seqlevels(nsegs) = as.character(1:22)
    kag = JaBbA:::karyograph(junctions = junc, tile = nsegs, label.edges = TRUE)
    expect_equal(length(kag$tile), 1144) ## removed garbage seqnames
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

#' xtYao #' Friday, Feb 26, 2021 10:39:37 AM
#' New testing data, rigma in 
jj = system.file("testing", "junctions.rds", package = "JaBbA")
cf = system.file("testing", "coverage.txt", package = "JaBbA")
ht = system.file("testing", "hets.txt", package = "JaBbA")
hr = fread(ht) %>% dt2gr

## default is boolean
jab = suppressWarnings(
    JaBbA(junctions = jj,
          coverage = cf,
          whitelist.junctions = whitelist.junctions,
          blacklist.coverage = blacklist.coverage,
          ## seg = segs,
          ## nseg = nsegs,
          ## strict = TRUE,
          slack.penalty = 10,
          hets = ht,
          tilim = 60,
          cfield = 'nudge',
          verbose = 2,
          outdir = './JaBbA.allin',
          overwrite = TRUE,
          ploidy=4.5,## preset HCC1954
          purity=1,
          epgap = 0.01,
          all.in = TRUE,
          ## juncs.uf = juncs.fn,
          tfield = 'nothing',
          nudge.balanced = TRUE,
          dyn.tuning = TRUE,
          max.na = 1)
)

## wj = readRDS(whitelist.junctions)
## mj = merge(jJ(juncs), jab$junctions)

## with iteration, linear penalty, no dynamic tuning
jab.reiterate = JaBbA(junctions = jj,
                      coverage = cf,
                      field = "ratio2",
                      blacklist.junctions = blacklist.junctions,
                      hets = hr,
                      slack.penalty = 10,
                      tilim = 60,
                      verbose = 2,
                      outdir = 'JaBbA.reiterate',
                      overwrite = TRUE,
                      reiterate=3,
                      ploidy=3.72,
                      purity=0.99,
                      loose.penalty.mode = 'linear',
                      epgap = 0.01,
                      dyn.tuning = FALSE,
                      max.na = 1)

## LP unit testing
## this test verfies consistency between LP and QP solutions
jab.lp = suppressWarnings(
    JaBbA(junctions = jj,
                    coverage = cf,
          whitelist.junctions = whitelist.junctions,
          blacklist.coverage = blacklist.coverage,
          ## seg = segs,
          ## nseg = nsegs,
          ## strict = TRUE,
          slack.penalty = 10,
          hets = ht,
          tilim = 60,
          cfield = 'nudge',
          verbose = 2,
          outdir = './JaBbA.lp',
          overwrite = TRUE,
          ploidy=4.57,## preset HCC1954
          purity=1,
          epgap = 1e-8
          all.in = TRUE,
          ## juncs.uf = juncs.fn,
          tfield = 'nothing',
          nudge.balanced = TRUE,
          dyn.tuning = TRUE,
          lp = TRUE,
          max.na = 1)
)

## for testing purposes, print out the exact output
print('jab cn')
print(list.expr(
    gr.string(sort(gr.stripstrand(jab$gr %Q% (strand=="+" & !is.na(cn)))), other.cols="cn")
))

print('jab junctions cn')
print(list.expr(values(jab$junctions$grl)$cn))

print('jab.reiterate cn')
print(list.expr(
    gr.string(sort(gr.stripstrand(jab.reiterate$gr %Q% (strand=="+" & !is.na(cn)))), other.cols="cn")
))

print('jab.reiterate junctions cn')
print(list.expr(values(jab.reiterate$junctions$grl)$cn))

print('jab.lp junctions cn')
print(list.expr(values(jab.lp$junctions$grl)$cn))

## print('jab.reiterate purity ploidy')
## print(paste(jab.reiterate$purity, jab.reiterate$ploidy))

## cn.cor.single = function(segs,
##                          cn.gs){
##     if (is.null(segs) | is.na(segs) | length(segs)==0){
##         return(as.numeric(NA))
##     }
##     bands.td = gTrack::karyogram()
##     bands.td$height=5
##     bands = bands.td@data
##     bands = grl.unlist(do.call(`GRangesList`, bands))
##     eligible = bands %Q% (stain != "acen") ## excluding CENTROMERE

##     ## reduce eligible region
##     rd.el = reduce(eligible + 1e4) - 1e4
##     rd.el.td = gTrack(rd.el)
    
##     cn.gs = cn.gs %*% rd.el ## select only overlaps

##     ov = gr2dt(gr.findoverlaps(cn.gs[, "cn"], segs[,"cn"]))
##     ov[, ":="(cn = segs$cn[subject.id],
##               gs.cn = cn.gs$cn[query.id],
##               gs.wd = width(cn.gs)[query.id])]
##     ov[!is.na(cn),
##        ":="(broken.into = .N,
##             inferred.cn = sum(cn*width)/sum(width)),
##        by="query.id"]

##     sp.cor = ov[!duplicated(query.id) & gs.cn<500, cor(inferred.cn, gs.cn, use="na.or.complete", method="spearman")]

##     return(sp.cor)
## }

## cn.gs = readRDS(system.file("extdata/jab.cn.gs.rds", package="JaBbA"))
## cn.gs.2 = readRDS(system.file("extdata/jab.cn.gs.2.rds", package="JaBbA"))
## cn.gs.reiterate = readRDS(system.file("extdata/jab.reiterate.cn.gs.rds", package="JaBbA"))
## cn.gs.reiterate.2 = readRDS(system.file("extdata/jab.reiterate.cn.gs.2.rds", package="JaBbA"))

test_that("JaBbA", {
    print("Comparing results from boolean mode without iteration:")
    expect_true(identical(jab$nodes$dt[!is.na(cn), cn], c(5, 3, 2, 4, 5, 3, 5, 3, 2, 3, 5)))
    ## expect_true((jab.cn.cor <<- pmax(
    ##                  cn.cor.single(jab$gr %Q% (strand=="+"), cn.gs),
    ##                  cn.cor.single(jab$gr %Q% (strand=="+"), cn.gs.2)
    ##              )) > 0.8,
    ##             info = print(jab.cn.cor))
    ## travis = c(3, 3, 3, 1, 3, 4, 2, 3, 2, 3, 2, 1, 2, 3, 3, 14, 11, 14, 14, 25, 29, 31, 2, 31, 31, 2, 31, 31, 2, 31, 31, 3, 31, 24, 7, 24, 29, 31, 2, 31, 33, 1, 33, 33, 16, 20, 30, 30, 33, 1, 33, 33, 1, 33, 32, 1, 32, 27, 6, 2, 25, 5, 4, 3, 1, 3, 4, 3, 3, 11, 4, 3, 0, 0)
    ## local = values(jab$junctions$grl)$cn
    ## cor(values(junc)$cool_cn, values(readRDS("JaBbA/junctions.rds"))$cn.jabba)
    expect_true(identical(jab$junctions$dt$cn, c(3, 2, 2, 1, 2, 4, 3, 2, 3, 3, 2, 2, 1, 2, 3)))
    ## expect_true(
    ##     identical(values(jab$junctions$grl)$cn,
    ##               c(3, 3, 3, 1, 3, 4, 2, 3, 2, 3, 2, 1, 2, 3, 3, 14, 11, 14, 14, 25, 28, 29, 31, 2, 31, 31, 2, 31, 31, 2, 31, 31, 3, 31, 24, 7, 24, 29, 31, 1, 31, 32, 1, 32, 32, 2, 32, 32, 16, 20, 29, 29, 33, 1, 33, 33, 1, 33, 32, 1, 32, 27, 6, 2, 25, 5, 4, 3, 1, 3, 4, 3, 3, 11, 3, 4, 0, 0)) |
    ##     identical(values(jab$junctions$grl)$cn,
    ##               c(3, 3, 3, 1, 3, 4, 2, 3, 2, 3, 2, 1, 2, 3, 3, 14, 11, 14, 14, 25, 29, 31, 2, 31, 31, 2, 31, 31, 2, 31, 31, 3, 31, 24, 7, 24, 29, 31, 2, 31, 33, 1, 33, 33, 17, 20, 29, 29, 33, 1, 33, 33, 1, 33, 32, 1, 32, 27, 6, 2, 25, 5, 4, 3, 1, 3, 4, 3, 3, 11, 4, 4, 0, 0)) |
    ##     identical(values(jab$junctions$grl)$cn,
    ##               c(3, 3, 3, 1, 3, 4, 2, 3, 2, 3, 2, 1, 2, 3, 3, 14, 11, 14, 14, 25, 29, 31, 2, 31, 31, 2, 31, 31, 2, 31, 31, 3, 31, 24, 7, 24, 29, 31, 2, 31, 33, 1, 33, 33, 16, 20, 30, 30, 33, 1, 33, 33, 1, 33, 32, 1, 32, 27, 6, 2, 25, 5, 4, 3, 1, 3, 4, 3, 3, 11, 4, 3, 0, 0)),
    ##     info = print(list.expr(values(jab$junctions$grl)$cn))
    ## )

    print("Comparing results from linear mode with iteration:")
    expect_true(identical(jab.reiterate$nodes$dt[!is.na(cn), cn], c(4, 3, 2, 3, 4, 3, 4, 3, 2, 3, 4)))
    ## expect_true((jab.reiterate.cn.cor <<- pmax(
    ##                  cn.cor.single(jab.reiterate$gr %Q% (strand=="+"), cn.gs.reiterate),
    ##                  cn.cor.single(jab.reiterate$gr %Q% (strand=="+"), cn.gs.reiterate.2)
    ##              )) > 0.8,
    ##             info = print(jab.reiterate.cn.cor))

    ## expect_true((jab.reiterate.cn.cor <<- cn.cor.single(jab.reiterate$gr %Q% (strand=="+"), cn.gs.reiterate)) > 0.8,
    ##             info = print(jab.reiterate.cn.cor))
    expect_true(identical(jab.reiterate$junctions$dt$cn, c(3, 1, 2, 1, 2, 3, 3, 1, 3, 3, 1, 2, 1, 2, 3)))
    ## expect_true(
    ##     identical(
    ##         values(jab.reiterate$junctions$grl)$cn,
    ##         c(2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 3, 2, 2, 3, 4, 4, 5, 5, 4, 3, 3, 2, 3, 2, 3, 3, 2, 2, 3, 3, 6, 9, 14, 17, 16, 16, 18, 19, 28, 29, 30, 29, 9, 29, 31, 31, 31, 32, 27, 6, 27, 31, 31, 31, 29, 27, 7, 24, 24, 24, 12, 6, 5, 4, 3, 3, 3, 3, 4, 4, 5, 7, 9, 7, 5, 4, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 9, 18, 31, 0, 0, 0, 0, 0)) |
    ##     identical(
    ##         values(jab.reiterate$junctions$grl)$cn,
    ##         c(2, 3, 3, 3, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 2, 2, 3, 4, 5, 4, 4, 5, 5, 4, 3, 2, 1, 2, 3, 3, 13, 11, 13, 13, 23, 28, 26, 6, 26, 31, 22, 19, 22, 31, 1, 31, 27, 5, 3, 24, 4, 3, 1, 3, 4, 3, 3, 3, 3, 10, 5, 10, 0, 0)),
    ##     info = print(list.expr(values(jab.reiterate$junctions$grl)$cn)))
    print("Comparing results from LP mode with L0 penalty")
    expect_true(identical(jab.lp$nodes$dt[!is.na(cn) & cn > 0, cn], c(5, 3, 2, 4, 5, 3, 5, 3, 2, 3, 5)))
    expect_true(identical(jab.lp$junctions$dt$cn, c(3, 2, 2, 1, 2, 4, 3, 2, 3, 3, 2, 2, 1, 2, 3)))
})
