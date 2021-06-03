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

## check that all junctions are incorporated even with garbage coverage

message("Testing incorporation of Tier 1 junctions even with ISM = TRUE")

dup.juncs.lp = suppressWarnings(
    jbaLP(gg.file = dup.juncs.gg,
         tilim = 60,
         tfield = "tier",
         cn.field = "cn",
         ism = TRUE,
         epgap = 1e-2,
         verbose = 2,
         return.type = "gGraph")
)

inp.juncs = readRDS(dup.juncs.grl)

## all should be incorporated
expect_equal(length(dup.juncs.lp$edges$dt[type == "ALT" & cn > 0, edge.id]),
             length(inp.juncs))

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
