library(JaBbA)
library(gUtils)
library(testthat)
jj = system.file("testing", "junctions.rds", package = "JaBbA")
sg = system.file("testing", "segs.rds", package = "JaBbA")
cf = system.file("testing", "coverage.txt", package = "JaBbA")

## save stuff as bw
cf.bw = system.file("testing", "coverage.bw", package = "JaBbA")
seg.bw = system.file("testing", "segs.bw", package = "JaBbA")

expected.cns = c(5,3,2,4,5,3,5,3,2,3,5)

test_that(desc = "Testing LP JaBbA with .bw coverage input",
          code = {
              jba.lp = suppressWarnings(JaBbA(
                  junctions = jj,
                  coverage = cf.bw,
                  seg = sg,
                  slack.penalty = 10,
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
                  max.na = 1))
              message(paste(jba.lp$nodes$dt$cn, collapse = "\n"))
              expect_equal(jba.lp$nodes$dt[cn > 0, cn],
                           expected.cns,
                           tolerance = 4)
          })

test_that(desc = "Testing LP JaBbA with .bw segment input",
          code = {
              jba.lp = suppressWarnings(JaBbA(
                  junctions = jj,
                  coverage = cf,
                  seg = seg.bw,
                  slack.penalty = 10,
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
                  max.na = 1))
              expect_equal(jba.lp$nodes$dt[cn > 0, cn],
                           expected.cns,
                           tolerance = 4)
          })

test_that(desc = "Testing LP JaBbA with .bw segment AND coverage input",
          code = {
              jba.lp = suppressWarnings(JaBbA(
                  junctions = jj,
                  coverage = cf.bw,
                  seg = seg.bw,
                  slack.penalty = 10,
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
                  max.na = 1))
              expect_equal(jba.lp$nodes$dt[cn > 0, cn],
                           expected.cns,
                           tolerance = 4)
          })
