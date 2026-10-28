d <- analysis_data(sim_trial(600, "S3", seed = 30))

test_that("the standard forest plot has all 13 subgroups", {
  res <- m_naive(d)
  expect_equal(nrow(res), 13)
  expect_true(all(is.finite(res$est)))
  expect_true(all(res$lo < res$est & res$est < res$hi))
})

test_that("interaction and biomarker tests give p-values", {
  p_int <- m_interaction(d)
  expect_named(p_int, sg_vars)
  expect_true(all(p_int >= 0 & p_int <= 1))
  p_bm <- m_biomarker(d)
  expect_named(p_bm, c("cutoff", "linear", "spline"))
  expect_true(all(p_bm >= 0 & p_bm <= 1))
})

test_that("the overall Cox model recovers the true HR in a large trial", {
  big <- analysis_data(sim_trial(6000, "S2", seed = 31))
  expect_equal(unname(overall_cox(big)["est"]), log(0.70), tolerance = 0.1)
})
