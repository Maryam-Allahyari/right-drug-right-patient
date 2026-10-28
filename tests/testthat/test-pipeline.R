# The same check as trial/scripts/01_validate_pipeline.R, on a small trial:
# simulate -> SDTM -> ADaM must give back exactly the OS and PFS the
# simulation defines, for every patient.

pipeline <- function(n, scenario, seed) {
  d <- sim_trial(n, scenario, seed = seed)
  sdtm <- write_sdtm(d)
  adsl <- derive_adsl(sdtm$dm, sdtm$ex, sdtm$ds, sdtm$sc, sdtm$rs, sdtm$ss)
  adrs <- derive_adrs(sdtm$rs, adsl)
  list(truth = pfs_os_truth(d) %>% mutate(USUBJID = sprintf("NSCLC01-%04d", id)),
       adsl = adsl, adtte = derive_adtte(adsl, adrs))
}

for (s in c("S1", "S4")) {
  test_that(paste("ADTTE matches the simulated OS and PFS in", s), {
    p <- pipeline(150, s, seed = 10 + match(s, c("S1", "S4")))
    cmp <- p$adtte %>%
      select(USUBJID, PARAMCD, AVAL, CNSR) %>%
      tidyr::pivot_wider(names_from = PARAMCD, values_from = c(AVAL, CNSR)) %>%
      inner_join(p$truth, by = "USUBJID")
    expect_equal(nrow(cmp), 150)
    expect_equal(cmp$AVAL_OS, cmp$OS_AVAL)
    expect_equal(cmp$CNSR_OS, 1 - cmp$os_event)   # CNSR is 1 for censored
    expect_equal(cmp$AVAL_PFS, cmp$PFS_AVAL)
    expect_equal(cmp$CNSR_PFS, 1 - cmp$pfs_event)
  })
}

test_that("ADSL has one row per patient, all in the ITT population", {
  p <- pipeline(100, "S2", seed = 20)
  expect_equal(nrow(p$adsl), 100)
  expect_false(anyDuplicated(p$adsl$USUBJID) > 0)
  expect_true(all(p$adsl$ITTFL == "Y"))
  expect_setequal(unique(p$adsl$REGION1), c("North America", "Europe", "Asia"))
})
