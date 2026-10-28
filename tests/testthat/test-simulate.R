test_that("the same seed gives the same trial", {
  expect_identical(sim_trial(200, "S2", seed = 1), sim_trial(200, "S2", seed = 1))
})

test_that("arms are exactly 1:1", {
  d <- sim_trial(301, "S1", seed = 2)
  expect_equal(nrow(d), 301)
  expect_lte(abs(sum(d$arm == 1) - sum(d$arm == 0)), 1)
})

test_that("true log HRs follow the scenario", {
  d <- sim_trial(500, "S1", seed = 3)
  expect_true(all(scenario_loghr("S1", d) == 0))
  expect_equal(unique(scenario_loghr("S2", d)), log(0.70))
  s3 <- scenario_loghr("S3", d)
  expect_equal(unique(s3[d$ecog == "0"]), log(0.55))
  expect_true(all(s3[d$ecog == "1"] == 0))
  # S4: HR 0.70 at the cut-off, smaller (better) for higher biomarker values
  s4 <- scenario_loghr("S4", d)
  expect_equal(s4[d$bmrk == 5], rep(log(0.70), sum(d$bmrk == 5)))
  expect_lt(cor(s4, d$bmrk), 0)
})

test_that("OS and PFS follow the plan's rules", {
  d <- sim_trial(400, "S4", seed = 4)
  tr <- pfs_os_truth(d)
  expect_true(all(tr$OS_AVAL >= 1 & tr$PFS_AVAL >= 1))
  # progression can't be after death, so PFS is never longer than OS
  expect_true(all(tr$PFS_AVAL <= tr$OS_AVAL))
  expect_true(all(tr$os_event %in% 0:1) && all(tr$pfs_event %in% 0:1))
  # patients alive are followed to the cut-off or drop-out, not just the last scan
  alive <- tr$os_event == 0
  expect_equal(tr$os_day[alive], tr$fu_end[alive])
})

test_that("the control arm is calibrated to about 11 months median OS", {
  d <- analysis_data(sim_trial(4000, "S1", seed = 5))
  km <- survfit(Surv(time, event) ~ 1, data = d[d$arm == 0, ])
  expect_equal(unname(summary(km)$table["median"]), 11, tolerance = 0.1)
})

test_that("S5 has no effect for 4 months and HR 0.60 after", {
  d <- analysis_data(sim_trial(20000, "S5", seed = 6))
  early <- transform(d, event = event * (time <= 4), time = pmin(time, 4))
  late <- subset(d, time > 4)
  expect_equal(unname(overall_cox(early)["est"]), 0, tolerance = 0.1)
  expect_equal(unname(overall_cox(late)["est"]), log(0.60), tolerance = 0.1)
})
