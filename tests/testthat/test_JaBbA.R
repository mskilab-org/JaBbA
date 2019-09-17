library(JaBbA)
library(gUtils)
library(gTrack)
library(testthat)

context('JaBbA')

print(sessionInfo())

juncs.fn = system.file("extdata", "junctions.vcf", package = 'JaBbA')
bedpe = system.file("extdata", "junctions.bedpe", package = 'JaBbA')
cov.fn = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')
segs = system.file("extdata", "segs.rds", package = 'JaBbA')

Sys.setenv(DEFAULT_GENOME = paste0(system.file(
               "extdata",
               package = 'JaBbA'),
               "/human_g1k_v37.chrom.sizes")) ## must do this           
print(Sys.getenv("DEFAULT_GENOME"))
print("This is the default seqlengths: ")
print(hg_seqlengths())
print("IS IT A VECTOR NOW??????!!!!!!")
print(is.vector(hg_seqlengths()))

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

junc = read.junctions(juncs.fn)
values(junc)$nudge = 0
junc = rep(junc, 2)

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
TILIM = 900
EPGAP = 0.95
nsegs = readRDS(segs)
nsegs$cn = 2

hets.gr = dt2gr(fread(hets))

test_that("karyograph", {
    kag = JaBbA:::karyograph(junctions = junc)
    expect_equal(length(kag$tile), 382)
    seqlevels(nsegs) = as.character(1:22)
    kag = JaBbA:::karyograph(junctions = junc, tile = nsegs, label.edges = TRUE)
    expect_equal(length(kag$tile), 1190)
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

## default is boolean
jab = JaBbA(junctions = junc,
            coverage = cov.fn,
            seg = segs,
            nseg = nsegs,
            strict = TRUE,
            slack.penalty = 1e4,
            hets = hets,
            tilim = TILIM,
            cfield = 'nudge',
            verbose = 2,
            outdir = 'JaBbA',
            overwrite = TRUE,
            ploidy=3.72,
            purity=NA,
            epgap = EPGAP,
            all.in = TRUE,
            juncs.uf = juncs.fn,
            tfield = 'nothing',
            nudge.balanced = TRUE,
            dyn.tuning = TRUE)

## with iteration, linear penalty, no dynamic tuning
jab.reiterate = JaBbA(junctions = juncs.fn,
                      coverage = cov.fn,
                      hets = hets.gr,
                      slack.penalty = 1e4,
                      tilim = TILIM,
                      verbose = 2,
                      outdir = 'JaBbA.reiterate',
                      overwrite = TRUE,
                      reiterate=3,
                      ploidy=3.72,
                      purity=0.99,
                      loose.penalty.mode = 'linear',
                      epgap = EPGAP,
                      dyn.tuning = FALSE)

print('jab cn')
print(class(jab))
print(jab$segstats)
print(list.expr(
    gr.string(sort(gr.stripstrand(jab$gr %Q% (strand=="+"))), other.cols="cn")
))

print('jab junctions cn')
print(list.expr(values(jab$junctions$grl)$cn))

## print('jab purity ploidy')
## print(paste(jab$purity, jab$ploidy))

print('jab.reiterate cn')
print(class(jab.reiterate))
print(jab.reiterate$segstats)
print(list.expr(
    gr.string(sort(gr.stripstrand(jab.reiterate$gr %Q% (strand=="+"))), other.cols="cn")
))

print('jab.reiterate junctions cn')
print(list.expr(values(jab.reiterate$junctions$grl)$cn))

## print('jab.reiterate purity ploidy')
## print(paste(jab.reiterate$purity, jab.reiterate$ploidy))

cn.cor.single = function(segs,
                         cn.gs){
    if (is.null(segs) || is.na(segs) || length(segs)==0){
        return(as.numeric(NA))
    }
    bands.td = gTrack::karyogram()
    bands.td$height=5
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    eligible = bands %Q% (stain != "acen") ## excluding CENTROMERE

    ## reduce eligible region
    rd.el = reduce(eligible + 1e4) - 1e4
    rd.el.td = gTrack(rd.el)
    
    cn.gs = cn.gs %*% rd.el ## select only overlaps

    ov = gr2dt(gr.findoverlaps(cn.gs[, "cn"], segs[,"cn"]))
    ov[, ":="(cn = segs$cn[subject.id],
              gs.cn = cn.gs$cn[query.id],
              gs.wd = width(cn.gs)[query.id])]
    ov[!is.na(cn),
       ":="(broken.into = .N,
            inferred.cn = sum(cn*width)/sum(width)),
       by="query.id"]

    sp.cor = ov[!duplicated(query.id) & gs.cn<500, cor(inferred.cn, gs.cn, use="na.or.complete", method="spearman")]

    return(sp.cor)
}

cn.gs = readRDS(system.file("extdata/jab.cn.gs.rds", package="JaBbA"))
cn.gs.2 = readRDS(system.file("extdata/jab.cn.gs.2.rds", package="JaBbA"))
cn.gs.reiterate = readRDS(system.file("extdata/jab.reiterate.cn.gs.rds", package="JaBbA"))
cn.gs.reiterate.2 = readRDS(system.file("extdata/jab.reiterate.cn.gs.2.rds", package="JaBbA"))

test_that("JaBbA", {
    print("Comparing results from boolean mode without iteration:")

    expect_true((jab.cn.cor <<- pmax(
                     cn.cor.single(jab$gr %Q% (strand=="+"), cn.gs),
                     cn.cor.single(jab$gr %Q% (strand=="+"), cn.gs.2)
                 )) > 0.8,
                info = print(jab.cn.cor))

    expect_true(identical(values(jab$junctions$grl)$cn,
                          c(3, 3, 3, 29, 3, 29, 28, 1, 28, 5, 28, 16, 12, 16, 16, 16, 16, 8, 16, 4, 3, 3, 0, 0)) |
                identical(values(jab$junctions$grl)$cn,
                          c(3, 3, 3, 29, 3, 29, 28, 1, 28, 5, 28, 16, 12, 16, 16, 16, 16, 8, 16, 4, 3, 3, 0, 0)),
                info = print(list.expr(values(jab$junctions$grl)$cn)))

    ## expect_true(identical(values(jab$junctions$grl)$cn, c(3, 3, 1, 6, 12, 17, 8, 3, 3, 29, 29, 28, 28, 28, 16, 16, 16, 16, 16, 4, 3, 3)) |
    ##             identical(values(jab$junctions$grl)$cn, c(3, 3, 3, 14, 20, 23, 28, 29, 30, 31, 28, 4, 28, 30, 31, 3, 31, 29, 2, 29, 4,
    ##                                                       29, 28, 7, 4, 3, 3, 4, 7, 5, 4, 3, 3, 2, 3, 3, 3)),
    ##             info = print(list.expr(values(jab$junctions$grl)$cn)))

    ## expect_true(abs(jab$ploidy - 3.50)<0.2,
    ##             info = print(jab$ploidy))

    ## expect_true(abs(jab$purity - .98)<0.02 |
    ##             abs(jab$purity -  1.000000)<0.01,
    ##             info = print(jab$purity))

    print("Comparing results from linear mode with iteration:")
    expect_true((jab.reiterate.cn.cor <<- pmax(
                     cn.cor.single(jab.reiterate$gr %Q% (strand=="+"), cn.gs.reiterate),
                     cn.cor.single(jab.reiterate$gr %Q% (strand=="+"), cn.gs.reiterate.2)
                 )) > 0.8,
                info = print(jab.reiterate.cn.cor))

    expect_true((jab.reiterate.cn.cor <<- cn.cor.single(jab.reiterate$gr %Q% (strand=="+"), cn.gs.reiterate)) > 0.8,
                info = print(jab.reiterate.cn.cor))

    expect_true(
        identical(
          values(jab.reiterate$junctions$grl)$cn,
            c(2, 3, 3, 3, 4, 3, 4, 3, 4, 4, 4, 4, 5, 3, 4, 2, 2, 4, 5, 4, 4, 4, 4, 3, 2, 1, 2, 3, 3, 12, 12, 12, 12, 21, 26, 27, 7, 27, 31, 1, 31, 32, 25, 13, 25, 31, 1, 31, 25, 7, 8, 17, 4, 3, 1, 3, 4, 3, 3, 9, 5, 6, 0, 0)) |
        identical(
            values(jab.reiterate$junctions$grl)$cn,
            c(2, 3, 3, 3, 4, 3, 4, 3, 4, 4, 4, 4, 5, 3, 5, 2, 2, 4, 5, 4, 4, 5, 4, 3, 2, 1, 2, 3, 21, 3, 3, 20, 16, 10, 16, 31, 2, 31, 33, 26, 12, 26, 23, 9, 6, 16, 4, 4, 3, 3, 3, 3, 17, 6)),
      info = print(list.expr(values(jab.reiterate$junctions$grl)$cn)))

    ## expect_true(identical(values(jab.reiterate$junctions$grl)$cn,
    ##                       c(4, 4, 1, 3, 10, 11, 5, 3, 6, 10, 19, 1, 1, 5, 1, 1, 2, 3, 3,
    ##                         3, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 2, 2, 3, 4, 5, 4, 4, 5, 5,
    ##                         4, 4, 4, 4, 4, 3, 2, 2, 3, 13, 13, 13, 23, 29, 26, 26, 31,
    ##                         22, 22, 32, 32, 31, 31, 27, 23, 23, 23, 4, 3, 3, 4, 3, 3)) |
    ##             identical(values(jab.reiterate$junctions$grl)$cn,
    ##                       c(2, 3, 3, 3, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 2, 2, 3, 4,
    ##                         5, 4, 4, 5, 5, 4, 4, 4, 4, 4, 3, 2, 1, 2, 3, 3, 13, 11, 13,
    ##                         13, 23, 29, 26, 6, 26, 31, 22, 19, 22, 32, 1, 32, 31, 1, 31,
    ##                         27, 5, 3, 23, 23, 1, 23, 4, 3, 1, 3, 4, 3, 3, 10, 5, 10, 0, 0)),
    ##             info = print(list.expr(values(jab.reiterate$junctions$grl)$cn)))

    ## expect_true(abs(jab.reiterate$ploidy - 3.62)<0.01 |
    ##             abs(jab.reiterate$ploidy - 3.51)<0.01, info = print(jab.reiterate$ploidy))

    ## expect_true(abs(jab.reiterate$purity - .99)<0.01 |
    ##             abs(jab.reiterate$purity - .99)<0.01,
    ##             info = print(jab.reiterate$purity))

})
