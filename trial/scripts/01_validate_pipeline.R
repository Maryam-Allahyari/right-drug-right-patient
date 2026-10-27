# Checks the SDTM -> ADaM pipeline against the simulated truth.
# For 20 trials per scenario (N = 300): simulate, write SDTM, derive
# ADSL/ADRS/ADTTE with admiral/admiralonco, then compare every patient's
# OS and PFS (days and censoring) with pfs_os_truth().
# If nothing differs, the big study can use analysis_data() directly
# instead of going through SDTM 4,000 times.
#
# Before running:
#   install.packages(c("admiral", "admiralonco", "dplyr", "tidyr"))
#
# From the repo root:
#   Rscript trial/scripts/01_validate_pipeline.R
#
# Output: trial/results/pipeline_vs_truth.csv (one row per trial) and a
# summary per scenario on screen. Takes a few minutes.

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

out <- "trial/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

check_trial <- function(scenario, rep, n = 300) {
  # seeds are kept apart from the ones the big study uses
  d <- sim_trial(n, scenario, seed = 1e5 + 100 * match(scenario, c("S1", "S2", "S3", "S4")) + rep)
  truth <- pfs_os_truth(d) %>% mutate(USUBJID = sprintf("NSCLC01-%04d", id))
  sdtm <- write_sdtm(d)

  adsl  <- derive_adsl(sdtm$dm, sdtm$ex, sdtm$ds, sdtm$sc, sdtm$rs, sdtm$ss)
  adrs  <- derive_adrs(sdtm$rs, adsl)
  adtte <- derive_adtte(adsl, adrs)

  cmp <- adtte %>%
    select(USUBJID, PARAMCD, AVAL, CNSR) %>%
    tidyr::pivot_wider(names_from = PARAMCD, values_from = c(AVAL, CNSR)) %>%
    inner_join(truth, by = "USUBJID")

  # CNSR is 1 for censored, so it should be the opposite of the event flag
  tibble(
    scenario = scenario, rep = rep, n = nrow(cmp),
    os_aval_diff  = sum(cmp$AVAL_OS != cmp$OS_AVAL),
    os_cnsr_diff  = sum(cmp$CNSR_OS != 1 - cmp$os_event),
    pfs_aval_diff = sum(cmp$AVAL_PFS != cmp$PFS_AVAL),
    pfs_cnsr_diff = sum(cmp$CNSR_PFS != 1 - cmp$pfs_event),
    os_events = sum(cmp$os_event), pfs_events = sum(cmp$pfs_event)
  )
}

# trials are independent, so they run in parallel
cl <- makeCluster(min(8, detectCores() - 2))
invisible(clusterCall(cl, function(p) .libPaths(p), .libPaths()))   # same packages as here
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(dplyr); library(admiral); library(admiralonco) })
  for (f in c("simulate", "write_sdtm", "derive_adsl", "derive_adrs", "derive_adtte"))
    source(file.path("trial/R", paste0(f, ".R")))
}))
clusterExport(cl, "check_trial")
jobs <- expand.grid(rep = 1:20, scenario = c("S1", "S2", "S3", "S4"), stringsAsFactors = FALSE)
res <- parLapply(cl, seq_len(nrow(jobs)), function(i, jobs) check_trial(jobs$scenario[i], jobs$rep[i]),
                 jobs = jobs) %>% bind_rows()
stopCluster(cl)

write.csv(res, file.path(out, "pipeline_vs_truth.csv"), row.names = FALSE)

res %>%
  group_by(scenario) %>%
  summarise(trials = n(), patients = sum(n),
            across(ends_with("_diff"), sum),
            mean_os_events = mean(os_events), mean_pfs_events = mean(pfs_events)) %>%
  print()
