library(JaBbA)
library(gUtils)
library(testthat)

context('JaBbA')
jj = system.file("testing", "junctions.rds", package = "JaBbA")
empty.jj = system.file("testing", "empty.junctions.rds", package = "JaBbA")
tier.jj = system.file("testing", "tier.juncs.rds", package = "JaBbA")
sg = system.file("testing", "segs.rds", package = "JaBbA")
cf = system.file("testing", "coverage.txt", package = "JaBbA")
empty.cf = system.file("testing", "empty.cov.txt", package = "JaBbA")
ht = system.file("testing", "hets.txt", package = "JaBbA")
new.ht = system.file("testing", "new.hets.txt", package = "JaBbA")
hr = fread(ht) %>% dt2gr

blacklist.junctions = system.file("extdata", "blacklist.junctions.rds", package = 'JaBbA')
whitelist.junctions = system.file("extdata", "whitelist.junctions.rds", package = 'JaBbA')
blacklist.coverage = system.file("extdata", "hg19.blacklist.coverage.rds", package = 'JaBbA')

expected.cns = c(5,3,2,4,5,3,5,3,2,3,5)

dup.juncs = readRDS(system.file("testing", "duplicate.juncs.rds", package = "JaBbA"))
dup.juncs.gg = system.file("testing", "duplicate.juncs.gg.rds", package = "JaBbA")
dup.juncs.grl = system.file("testing", "duplicate.juncs.grl.rds", package = "JaBbA")

message("Testing duplicate breakpoint detection")
vanilla.op = detect_duplicate_breakpoints(dup.juncs, tfield = "tier")

expect_setequal(vanilla.op$dt$edge.id, c(306, 374)) ## test that correct edges are found

message("Testing duplicate breakpoints without tier field")
no.tier = detect_duplicate_breakpoints(dup.juncs, tfield = "dummy")

expect_equal(length(no.tier), 0)

message("Testing duplicate breakpoint detection with only REF edges")
ref.only = detect_duplicate_breakpoints(dup.juncs[type == "REF"], tfield = "tier")

expect_equal(length(ref.only), 0)

## check that fix.thres works as expected for jbaLP
kag.sg = readRDS(system.file("testing", "fix.thres.kag.rds", package = "JaBbA"))

test_that(desc = "test that treemem works as expected",
          code = {
              expect_error(
                  object = { JaBbA:::jbaLP(gg = kag.sg,
                                          cn.field = "cnmle",
                                          var.field = "loess.var",
                                          bins.field = "nbins",
                                          epgap = 1e-6,
                                          tilim = 100,
                                          fix.thres = 10,
                                          lambda = 100,
                                          max.mem = 1)
                  },
                  regexp = "Not enough memory")
          })

test_that(desc = "fit.thres in jbaLP",
          code = {
              res = JaBbA:::jbaLP(gg = kag.sg,
                    cn.field = "cnmle",
                    var.field = "loess.var",
                    bins.field = "nbins",
                    epgap = 1e-6,
                    tilim = 100,
                    fix.thres = 10,
                    lambda = 100)
              res.dt = as.data.table(res$segstats)[fixed == TRUE,]

              ## test that cn is fixed to cnmle
              expect_equal(res.dt$cn, res.dt$cn.old)

              ## test that weights are greater than fix.thres
              expect_true(all(res.dt$unfixed.weight > 1e3))
          })

test_that(desc = "minimum fit.thres warning",
          code = {
              expect_warning(
              object = {
                  JaBbA:::jbaLP(gg = kag.sg,
                    cn.field = "cnmle",
                    var.field = "loess.var",
                    bins.field = "nbins",
                    epgap = 1e-6,
                    tilim = 100,
                    fix.thres = 1,
                    lambda = 100)
              },
              regexp = "Small value for fix.thres selected")
          })

test_that(desc = "no beta in karyograph warning",
          code = {
              expect_warning(
              object = {
                  kag.nobeta = kag.sg$copy
                  kag.nobeta$set(beta = NA)
                  JaBbA:::jbaLP(gg = kag.nobeta,
                    cn.field = "cnmle",
                    var.field = "loess.var",
                    bins.field = "nbins",
                    epgap = 1e-6,
                    tilim = 100,
                    fix.thres = 100,
                    lambda = 100)
              })
          })


test_that(desc = "fit.thres default setting",
          code = {
              res = JaBbA:::jbaLP(gg = kag.sg,
                                  cn.field = "cnmle",
                                  var.field = "loess.var",
                                  bins.field = "nbins",
                                  epgap = 1e-6,
                                  tilim = 100,
                                  fix.thres = -1,
                                  lambda = 100)
              expect_false("fixed" %in% names(values(res$segstats)))
              expect_false("unfixed.weight" %in% names(values(res$segstats)))
              expect_false("unfixed.cn" %in% names(values(res$segstats)))
          })



## check that all junctions are incorporated even with garbage coverage

test_that(desc = "Test incorporation of Tier 1 junctions even with ISM = TRUE",
          code = {

              dup.juncs.lp = suppressWarnings(
                  jbaLP(gg.file = dup.juncs.gg,
                        tilim = 60,
                        tfield = "tier",
                        cn.field = "cn",
                        ism = TRUE,
                        epgap = 1e-6,
                        verbose = 2,
                        return.type = "gGraph")
              )

              inp.juncs = readRDS(dup.juncs.grl)

              ## all should be incorporated
              expect_equal(length(dup.juncs.lp$edges$dt[type == "ALT" & cn > 0, edge.id]),
                           length(inp.juncs))

          })

test_that(desc = "Testing LP JaBbA without input junctions",
          code = {

              expect_warning(JaBbA(
                  junctions = "",
                  coverage = cf,
                  slack.penalty = 100,
                  tilim = 60,
                  cfield = 'ratio',
                  verbose = 2,
                  outdir = 'JaBbA.lp',
                  overwrite = TRUE,
                  ploidy=4.5,## preset HCC1954
                  purity=1,
                  epgap = 0.01,
                  all.in = TRUE,
                  tfield = 'nothing',
                  nudge.balanced = TRUE,
                  dyn.tuning = TRUE,
                  lp = TRUE,
                  ism = FALSE,
                  max.na = 1),
                  regexp = "no junction file is given")
          })

test_that(desc = "Testing LP JaBbA with wrong coverage field",
          code = {

              expect_warning(JaBbA(
                  junctions = "",
                  coverage = cf,
                  slack.penalty = 100,
                  tilim = 60,
                  cfield = 'ratio',
                  field = "bad.coverage.field",
                  verbose = 2,
                  outdir = 'JaBbA.lp',
                  overwrite = TRUE,
                  ploidy=4.5,## preset HCC1954
                  purity=1,
                  epgap = 1e-6,
                  all.in = TRUE,
                  tfield = 'nothing',
                  nudge.balanced = TRUE,
                  dyn.tuning = TRUE,
                  lp = TRUE,
                  ism = FALSE,
                  max.na = 1),
                  regexp = "bad.coverage.field not found in coverage GRanges metadata")
          })

test_that(desc = "Testing vanilla JaBbA LP",
          code = {

              jab.lp = suppressWarnings(
                  JaBbA(junctions = jj,
                        coverage = cf,
                        whitelist.junctions = whitelist.junctions,
                        blacklist.coverage = blacklist.coverage,
                        slack.penalty = 100,
                        hets = ht,
                        tilim = 60,
                        cfield = 'nudge',
                        verbose = 2,
                        outdir = 'JaBbA.lp',
                        overwrite = TRUE,
                        ploidy=4.57,## preset HCC1954
                        purity=1,
                        epgap = 1e-6,
                        all.in = TRUE,
                        tfield = 'nothing',
                        nudge.balanced = TRUE,
                        dyn.tuning = TRUE,
                        lp = TRUE,
                        ism = FALSE,
                        max.na = 1)
              )
              expect_equal(jab.lp$nodes$dt[cn > 0, cn], expected.cns, tolerance = 1.1)
          })

test_that(desc = "Testing vanilla JaBbA LP with gurobi",
          code = {
              if (requireNamespace("gurobi", quietly = TRUE)) {
                  jab.lp = suppressWarnings(
                      JaBbA(junctions = jj,
                            coverage = cf,
                            whitelist.junctions = whitelist.junctions,
                            blacklist.coverage = blacklist.coverage,
                            slack.penalty = 100,
                            hets = ht,
                            tilim = 60,
                            cfield = 'nudge',
                            verbose = 2,
                            outdir = 'JaBbA.lp',
                            overwrite = TRUE,
                            ploidy=4.57,## preset HCC1954
                            purity=1,
                            epgap = 1e-6,
                            all.in = TRUE,
                            tfield = 'nothing',
                            nudge.balanced = TRUE,
                            use.gurobi = TRUE,
                            dyn.tuning = FALSE,
                            lp = TRUE,
                            ism = FALSE,
                            max.na = 1)
                  )
                  ## only check for chromosome 12...
                  expect_equal(jab.lp$nodes$dt[seqnames == "12", cn], expected.cns, tolerance = 0.5)
              } else {
                  expect_error(object = {
                      jab.lp = suppressWarnings(
                          JaBbA(junctions = jj,
                                coverage = cf,
                                whitelist.junctions = whitelist.junctions,
                                blacklist.coverage = blacklist.coverage,
                                slack.penalty = 100,
                                hets = ht,
                                tilim = 60,
                                cfield = 'nudge',
                                verbose = 2,
                                outdir = 'JaBbA.lp',
                                overwrite = TRUE,
                                ploidy=4.57,## preset HCC1954
                                purity=1,
                                epgap = 1e-6,
                                all.in = TRUE,
                                tfield = 'nothing',
                                nudge.balanced = TRUE,
                                use.gurobi = TRUE,
                                dyn.tuning = TRUE,
                                lp = TRUE,
                                ism = FALSE,
                                max.na = 1)
                      )
                  })
              }
          })


test_that(desc = "Testing JaBbA LP with ISM",
          code = {
              jab.ism = suppressWarnings(
                  JaBbA(junctions = jj,
                        coverage = cf,
                        whitelist.junctions = whitelist.junctions,
                        blacklist.coverage = blacklist.coverage,
                        slack.penalty = 100,
                        hets = ht,
                        tilim = 60,
                        cfield = 'nudge',
                        verbose = 2,
                        outdir = 'JaBbA.mem',
                        overwrite = TRUE,
                        ploidy=4.57,## preset HCC1954
                        purity=1,
                        epgap = 1e-6,
                        all.in = TRUE,
                        tfield = 'nothing',
                        nudge.balanced = TRUE,
                        dyn.tuning = TRUE,
                        lp = TRUE,
                        ism = TRUE,
                        max.na = 1)
              )

              expect_equal(jab.ism$nodes$dt[cn > 0, cn], expected.cns, tolerance = 1.1)
          })

test_that("Testing JaBbA LP with tier 1 junctions and uniformly random coverage",
          code = {
              jab.tier = suppressWarnings(
                  JaBbA(junctions = tier.jj,
                        coverage = empty.cf,
                        seg = sg,
                        whitelist.junctions = whitelist.junctions,
                        blacklist.coverage = blacklist.coverage,
                        slack.penalty = 10,
                        hets = ht,
                        tilim = 60,
                        field = 'empty.ratio', ## random uniform
                        cfield = 'nudge',
                        verbose = 2,
                        outdir = 'JaBbA.tier',
                        overwrite = TRUE,
                        ploidy=4.57,## preset HCC1954
                        purity=1,
                        epgap = 1e-6,
                        all.in = TRUE,
                        tfield = 'tier',
                        nudge.balanced = TRUE,
                        dyn.tuning = TRUE,
                        lp = TRUE,
                        ism = TRUE,
                        max.na = 1)
              )

              ## expect five ALT junctions
              expect_equal(jab.tier$edges$dt[type == "ALT" & cn > 0, .N], 5)
          })

test_that("Testing JaBbA LP with reiteration and empty junctions file",
          code = {
              jab.empty = suppressWarnings(
                  JaBbA(junctions = empty.jj,
                        juncs.uf = jj,
                        reiterate = 1,
                        coverage = cf,
                        whitelist.junctions = whitelist.junctions,
                        blacklist.coverage = blacklist.coverage,
                        slack.penalty = 10,
                        rescue.window = 1e4,
                        hets = ht,
                        tilim = 60,
                        cfield = 'nudge',
                        verbose = 2,
                        outdir = 'JaBbA.empty.jj',
                        overwrite = TRUE,
                        ploidy=4.57,## preset HCC1954
                        purity=1,
                        epgap = 1e-6,
                        all.in = TRUE,
                        tfield = 'nothing',
                        nudge.balanced = TRUE,
                        dyn.tuning = TRUE,
                        lp = TRUE,
                        ism = FALSE,
                        max.na = 1)
              )

              expect_equal(jab.empty$nodes$dt[cn > 0, cn], expected.cns, tolerance = 1.1)
          })

test_that("Testing JaBbA LP with empty junctions file and tiers in juncs.uf",
          code = {
              jab.empty.2 = suppressWarnings(
                  JaBbA(junctions = empty.jj,
                        juncs.uf = tier.jj,
                        reiterate = 1,
                        coverage = cf,
                        whitelist.junctions = whitelist.junctions,
                        blacklist.coverage = blacklist.coverage,
                        slack.penalty = 10,
                        rescue.window = 1e4,
                        hets = ht,
                        tilim = 60,
                        cfield = 'nudge',
                        verbose = 2,
                        outdir = 'JaBbA.empty.jj',
                        overwrite = TRUE,
                        ploidy=4.57,## preset HCC1954
                        purity=1,
                        epgap = 1e-6,
                        all.in = TRUE,
                        tfield = 'tier',
                        nudge.balanced = TRUE,
                        dyn.tuning = TRUE,
                        lp = TRUE,
                        ism = FALSE,
                        max.na = 1)
              )
              ## make sure that ALT edges are actually incorporated
              expect_true(jab.empty.2$edges$dt[type == "ALT" & cn > 0, .N] > 0)
          })


