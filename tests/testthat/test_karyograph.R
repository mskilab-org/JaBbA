library(testthat)

context('JaBbA')

## for testing any function upstream of junction balance
invalid.jjs = system.file("extdata", "invalid.junctions.rds", package = "JaBbA")

test_that(desc = "test filtering of junctions with invalid endpoints",
          code = {
              jj = readRDS(invalid.jjs)
              filtered.jj = JaBbA:::filter_oob_junctions(jj)
              expect_true(length(jj) > length(filtered.jj))
              expect_true(all(start(stack(filtered.jj)) >= 1))
              expect_warning(
                  object = {filtered.jj = JaBbA:::filter_oob_junctions(jj)},
                  regexp = "Some junction breakpoints are < 1"
              )})
          
        
kag = system.file("extdata", "fix.thres.kag.rds", package = "JaBbA")

test_that(desc = "test ppgrid with various perturbations to input",
          code = {
              kag.gg = readRDS(kag)
              vanilla.segstats = kag.gg$gr
              expect_true({
                  res = ppgrid(vanilla.segstats,
                               purity.min = 1,
                               ploidy.min = 1,
                               purity.max = 1,
                               ploidy.max = 5,
                               plot = F,
                               verbose = F)
                  nrow(res) >= 1})

              expect_true({
                  res = ppgrid(vanilla.segstats,
                               purity.min = NA,
                               ploidy.min = NA,
                               purity.max = NA,
                               ploidy.max = NA,
                               plot = F,
                               verbose = F)
                  nrow(res) >= 1})
              
              nomean.segstats = kag.gg$gr
              nomean.segstats$mean = NULL
              expect_error(ppgrid(nomean.segstats,
                                  purity.min = 1,
                                  ploidy.min = 0.1,
                                  purity.max = 1,
                                  ploidy.max = 3,
                                  plot = F,
                                  verbose = F))
          })
