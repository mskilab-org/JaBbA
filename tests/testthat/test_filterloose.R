library(JaBbA)
library(gUtils)
library(gGnome)
library(testthat)

jj = system.file("extdata", "junctions.rds", package = "JaBbA")
sg = system.file("extdata", "segs.rds", package = "JaBbA")
cf = system.file("extdata", "coverage.txt", package = "JaBbA")

test_that(desc = "Testing LP JaBbA with .bw coverage input",
          code = {
              expect_message({
                  jba.lp = suppressWarnings(JaBbA(
                      junctions = jj,
                      coverage = cf,
                      seg = sg,
                      slack.penalty = 10,
                      tilim = 60,
                      cfield = 'ratio',
                      reiterate = 2,
                      verbose = 2,
                      outdir = 'JaBbA.filter',
                      overwrite = TRUE,
                      ploidy=4.5,
                      purity=1,
                      epgap = 0.01,
                      all.in = TRUE,
                      tfield = NULL,
                      nudge.balanced = TRUE,
                      dyn.tuning = TRUE,
                      lp = TRUE,
                      ism = FALSE,
                      max.na = 1,
                      filter_loose = TRUE))
                  },
                  regexp = "Starting loose end annotation")
          })

test_that(desc = "Testing LP JaBbA with .bw coverage input",
          code = {
              expect_message({
                  jba.lp = suppressWarnings(JaBbA(
                      junctions = jj,
                      coverage = cf,
                      seg = sg,
                      slack.penalty = 10,
                      tilim = 60,
                      cfield = 'ratio',
                      verbose = 2,
                      outdir = 'JaBbA.rescue',
                      overwrite = TRUE,
                      ploidy=4.5,
                      purity=1,
                      epgap = 0.01,
                      all.in = TRUE,
                      tfield = NULL,
                      nudge.balanced = TRUE,
                      dyn.tuning = TRUE,
                      lp = TRUE,
                      ism = FALSE,
                      max.na = 1,
                      rescue.all = FALSE,
                      filter_loose = FALSE))
                  },
                  regexp = "Skipping loose end annotation")
          })

test_that(desc = "Testing LP JaBbA with .bw coverage input",
          code = {
              expect_message({
                  jba.lp = suppressWarnings(JaBbA(
                      junctions = jj,
                      coverage = cf,
                      seg = sg,
                      slack.penalty = 10,
                      tilim = 60,
                      cfield = 'ratio',
                      verbose = 2,
                      outdir = 'JaBbA.nofilter',
                      overwrite = TRUE,
                      ploidy=4.5,
                      purity=1,
                      epgap = 0.01,
                      all.in = TRUE,
                      tfield = NULL,
                      nudge.balanced = TRUE,
                      dyn.tuning = TRUE,
                      lp = TRUE,
                      ism = FALSE,
                      max.na = 1))
                  },
                  regexp = "Skipping loose end annotation")
          })
              
