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


message("Testing vanilla JaBbA LP")
jab.lp = suppressWarnings(
    JaBbA(junctions = jj,
          coverage = cf,
          whitelist.junctions = whitelist.junctions,
          blacklist.coverage = blacklist.coverage,
          slack.penalty = 10,
          hets = ht,
          tilim = 60,
          cfield = 'nudge',
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
          ism = FALSE)
)

expect_equal(jab.lp$nodes$dt[cn > 0, cn], expected.cns)

message("Testing set max.mem parameter")
jab.mem = suppressWarnings(
    JaBbA(junctions = jj,
          coverage = cf,
          whitelist.junctions = whitelist.junctions,
          blacklist.coverage = blacklist.coverage,
          slack.penalty = 10,
          hets = ht,
          tilim = 60,
          cfield = 'nudge',
          verbose = 2,
          outdir = 'JaBbA.mem',
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
          max.mem = 4)
)

expect_equal(jab.mem$nodes$dt[cn > 0, cn], expected.cns)

message("Testing JaBbA LP with ISM")
jab.ism = suppressWarnings(
    JaBbA(junctions = jj,
          coverage = cf,
          whitelist.junctions = whitelist.junctions,
          blacklist.coverage = blacklist.coverage,
          slack.penalty = 10,
          hets = ht,
          tilim = 60,
          cfield = 'nudge',
          verbose = 2,
          outdir = 'JaBbA.mem',
          overwrite = TRUE,
          ploidy=4.5,## preset HCC1954
          purity=1,
          epgap = 0.01,
          all.in = TRUE,
          tfield = 'nothing',
          nudge.balanced = TRUE,
          dyn.tuning = TRUE,
          lp = TRUE,
          ism = TRUE)
)

expect_equal(jab.ism$nodes$dt[cn > 0, cn], expected.cns)

message("Testing JaBbA LP with tier 1 junctions and uniformly random coverage")
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
          ploidy=4.5,## preset HCC1954
          purity=1,
          epgap = 0.01,
          all.in = TRUE,
          tfield = 'tier',
          nudge.balanced = TRUE,
          dyn.tuning = TRUE,
          lp = TRUE,
          ism = TRUE)
)

## expect five ALT junctions
expect_equal(jab.tier$edges$dt[type == "ALT" & cn > 0, .N], 5)

message("Testing JaBbA LP with reiteration and empty junctions file")
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
          ploidy=4.5,## preset HCC1954
          purity=1,
          epgap = 0.01,
          all.in = TRUE,
          tfield = 'nothing',
          nudge.balanced = TRUE,
          dyn.tuning = TRUE,
          lp = TRUE,
          ism = FALSE)
)

expect_equal(jab.empty$nodes$dt[cn > 0, cn], expected.cns)

message("Testing JaBbA LP with empty junctions file and tiers in juncs.uf")
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
          ploidy=4.5,## preset HCC1954
          purity=1,
          epgap = 0.01,
          all.in = TRUE,
          tfield = 'tier',
          nudge.balanced = TRUE,
          dyn.tuning = TRUE,
          lp = TRUE,
          ism = FALSE)
)

expect_equal(jab.empty.2$nodes$dt[cn > 0, cn], expected.cns)
