library(JaBbA)
library(gUtils)
library(gTrack)
library(testthat)

context('JaBbA')

juncs.fn = system.file("extdata", "junctions.vcf", package = 'JaBbA')
bedpe = system.file("extdata", "junctions.bedpe", package = 'JaBbA')
cov.fn = system.file("extdata", "coverage.txt", package = 'JaBbA')
hets = system.file("extdata", "hets.txt", package = 'JaBbA')
segs = system.file("extdata", "segs.rds", package = 'JaBbA')

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
    expect_equal(ncol(values(ram)), 30)
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

## default is boolean
jab = JaBbA(junctions = junc,
            coverage = cov.fn,
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
                      overwrite = TRUE,
                      reiterate=3,
                      ploidy=3.72,
                      purity=0.99,
                      loose.penalty.mode = 'linear',
                      epgap = EPGAP,
                      dyn.tuning = FALSE)

print('jab cn')
print(list.expr(jab$segstats$cn))

print('jab junctions cn')
print(list.expr(values(jab$junctions)$cn))

print('jab purity ploidy')
print(paste(jab$purity, jab$ploidy))

print('jab.reiterate cn')
print(list.expr(jab.reiterate$segstats$cn))

print('jab.reiterate junctions cn')
print(list.expr(values(jab.reiterate$junctions)$cn))

print('jab.reiterate purity ploidy')
print(paste(jab.reiterate$purity, jab.reiterate$ploidy))

cn.cor.single = function(segs,
                         cn.gs,
                         plot = FALSE){
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

    ## tru.pl = sum(cn.gs$cn * (width(cn.gs)/1e6)) / sum(width(cn.gs)/1e6)
    ## inferred.pl = sum(segs$cn * (width(segs)/1e6)) / sum(width(segs)/1e6)

    ov = gr2dt(gr.findoverlaps(cn.gs[, "cn"], segs[,"cn"]))
    ov[, ":="(cn = segs$cn[subject.id],
              gs.cn = cn.gs$cn[query.id],
              gs.wd = width(cn.gs)[query.id])]
    ov[!is.na(cn),
       ":="(broken.into = .N,
            inferred.cn = sum(cn*width)/sum(width)),
       by="query.id"]

    ## fold change compared to genome avg
    ## ov[, ":="(inferred.fc = inferred.cn/inferred.pl,
    ##           gs.fc = gs.cn/tru.pl)]

    sp.cor = ov[!duplicated(query.id) & gs.cn<500, cor(inferred.cn, gs.cn, use="na.or.complete", method="spearman")]
    
    if (plot){
        cov.td = gTrack(cov, y.field="ratio", circles=TRUE, lwd.border=0.3)
        segs.td = gTrack(segs, y.field="cn", name="inferred")
        gs.td = gTrack(cn.gs, y.field="cn", name="gs")
        
        ppdf(plot(c(cov.td, gs.td, segs.td), "19"))

        ppdf(plot(c(cov.td,
                    gs.td,
                    segs.td,
                    gTrack(win, height=5, col="red"),
                    gTrack(dt2gr(ov), col="green", name="overlaps", y.field="inferred.cn")),
                  win+1e4))

        ov[abs(inferred.cn-gs.cn)>3, table(gs.wd<1e4)]

        kag = readRDS("~/projects/CellLines/Flow.redo/JaBbA.boolean/JaBbA/HCC1143_nygc/debug.seqz.pl/karyograph.rds")
        bad.regions = kag$segstats %Q% (bad==TRUE)
        bad.regions$col="purple"
        bad.td = gTrack(unname(bad.regions), y.field="nbins.nafrac", name="bad", height=5)

        ov[abs(inferred.cn-gs.cn)>3 & gs.wd>=1e4]
        ov[abs(inferred.fc-gs.fc)>3 & gs.wd>=1e4]
        ov[abs(inferred.cn-gs.cn)>3 & gs.wd<1e4]
        
        win2 = dt2gr(ov[abs(inferred.fc-gs.fc)>3 & gs.wd>=1e4])
        win2 = dt2gr(ov[abs(inferred.cn-gs.cn)>3 & gs.wd<1e4])
        win2 = dt2gr(ov[abs(inferred.cn-gs.cn)>1 & gs.wd>=1e4])
        hood = new.gg$hood(win2, 1e5)

        ppdf(plot(
            c(bands.td,
              cov.td,
              bad.td,
              gs.td,
              ## segs.td,
              gg$td,
              new.gg$td,
              gTrack(win2, height=5, col="red"),
              gTrack(dt2gr(ov), col="green", name="overlaps", y.field="inferred.cn", height=5)),
            streduce(win2, 5e3)[11:20]), width=16)
        
        require(ggplot2)

        ppdf(print(
            ov %>%
            ## subset(!duplicated(query.id) & gs.cn<500 & gs.wd<1e4) %>%
            subset(!duplicated(query.id) & gs.cn<500 & gs.wd>=1e4) %>%
            ggplot() +
            geom_point(aes(x = inferred.cn,
                           y = gs.cn),
                       shape=1) +
            geom_abline(slope=1, intercept=0,
                        color="salmon", linetype="dashed")
            ## geom_point(aes(x = inferred.fc,
            ##                y = gs.fc),
            ##            shape=2)
        ))
        
    }

    return(sp.cor)
}

cn.gs = readRDS(system.file("extdata/jab.cn.gs.rds", package="JaBbA"))
cn.gs.reiterate = readRDS(system.file("extdata/jab.reiterate.cn.gs.rds", package="JaBbA"))

test_that("JaBbA", {
    print("Comparing results from boolean mode without iteration:")
    expect_true(cn.cor.single(jab$segstats, cn.gs)>0.9)
    ## expect_true( ## accounting for difference in run speed local vs travis
    ##     identical(jab$segstats$cn,
    ##               c(4, 3, 3, 1, 3, 29, 32, 29, 28, 34, 28, 16, 33, 16, 24, 16, 4, 3, 4, 4, 3, 3, 1, 3, 29, 32, 29, 28, 34, 28, 16, 33, 16, 24, 16, 4, 3, 4, 1, 26, 1, 26, 1, 1)) |
    ##     identical(jab$segstats$cn,
    ##               c(3, 3, 1, 3, 10, 15, 20, 24, 28, 29, 30, 31, 32, 29, 30, 31, 32, 34, 32, 39, 28, 33, 28, 27, 25, 16, 5, 4, 3, 4, 8, 9, 6, 4, 3, 4, 3, 4, 3, 4, 3, 3, 1, 3, 10, 15, 20, 24, 28, 29, 30, 31, 32, 29, 30, 31, 32, 34, 32, 39, 28, 33, 28, 27, 25, 16, 5, 4, 3, 4, 8, 9, 6, 4, 3, 4, 3, 4, 3, 4, 1, 1, 2, 9, 1, 1, 2, 1, 1, 1, 7, 5, 5, 4, 4, 1, 1, 1, 1, 1, 1, 1, 7, 1, 1, 4, 1, 1, 1, 1, 7, 5, 5, 4, 4, 1, 1, 1, 1, 1, 1, 1, 7, 1, 1, 4, 1, 1, 1, 1, 1, 1, 2, 9, 1, 1, 2, 1, 1, 1)),
    ##     info = print(list.expr(jab$segstats$cn)))

    expect_true(identical(values(jab$junctions)$cn, c(2, 12, 3, 6, 17, 8, 1)) |
                identical(values(jab$junctions)$cn, c(2, 2, 5, 3, 11)) |
                identical(values(jab$junctions)$cn, c(2, 2, 4, 3, 11)),
                info = print(list.expr(values(jab$junctions)$cn)))

    expect_true(abs(jab$ploidy - 3.60)<0.1 |
                abs(jab$ploidy - 3.50)<0.1,
                info = print(jab$ploidy))

    expect_true(abs(jab$purity - .98)<0.01 |
                abs(jab$purity -  1.000000)<0.01,
                info = print(jab$purity))

    print("Comparing results from linear mode with iteration:")
    expect_true(cn.cor.single(jab.reiterate$segstats, cn.gs.reiterate)>0.9)
    ## expect_true(identical(jab.reiterate$segstats$cn,
    ##                       c(4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 34, 33, 27, 31, 32, 31, 41, 23, 32, 33, 32, 31, 32, 27, 22, 4, 3, 4, 10, 3, 4, 4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 34, 33, 27, 31, 32, 31, 41, 23, 32, 33, 32, 31, 32, 27, 22, 4, 3, 4, 10, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)) |
    ##             identical(jab.reiterate$segstats$cn,
    ##                       c(4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 27, 31, 40, 23, 32, 31, 32, 27, 21, 4, 3, 4, 9, 3, 4, 3, 4, 2, 4, 3, 2, 3, 3, 1, 3, 13, 24, 13, 23, 27, 33, 27, 31, 40, 23, 32, 31, 32, 27, 21, 4, 3, 4, 9, 3, 1, 1)),
    ##             info = print(list.expr(jab.reiterate$segstats$cn)))

    expect_true(identical(values(jab.reiterate$junctions)$cn,
                          c(1, 2, 1, 2, 10, 11, 4, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 6, 0, 0, 0, 0, 0, 0, 0, 1, 9, 0, 18, 0, 0, 0, 1, 0, 0, 0, 1, 5, 1)) |
                identical(values(jab.reiterate$junctions)$cn,
                          c(1, 2, 1, 2, 10, 11, 4, 6, 6, 0, 9, 17, 1, 5, 1)),
                info = print(list.expr(values(jab.reiterate$junctions)$cn)))

    expect_true(abs(jab.reiterate$ploidy - 3.62)<0.01 |
                abs(jab.reiterate$ploidy - 3.51)<0.01, info = print(jab.reiterate$ploidy))

    expect_true(abs(jab.reiterate$purity - .99)<0.01 |
                abs(jab.reiterate$purity - .99)<0.01,
                info = print(jab.reiterate$purity))
})
