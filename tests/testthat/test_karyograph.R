library(testthat)

context('JaBbA')

## for testing any function upstream of junction balance
invalid.jjs = system.file("testing", "invalid.junctions.rds", package = "JaBbA")

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
          
        

