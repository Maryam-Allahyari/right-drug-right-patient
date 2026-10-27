# The true subgroup HRs for each scenario, from one huge simulated trial
# (200,000 patients). With that many patients the estimates are practically
# exact, and they use the same Cox model as the forest plot, so every
# method is compared with the right target.
#
# From the repo root:
#   Rscript trial/scripts/02_truth.R
#
# Output: trial/outputs/truth.rds and a table of true HRs on screen.

suppressPackageStartupMessages(library(dplyr))
source("trial/R/simulate.R")
source("trial/R/methods.R")

dir.create("trial/outputs", recursive = TRUE, showWarnings = FALSE)

truth <- lapply(c("S1", "S2", "S3", "S4"), function(s) {
  d <- analysis_data(sim_trial(2e5, s, seed = 777))
  bind_rows(
    tibble(subgroup = "overall", est = overall_cox(d)[["est"]]),
    m_naive(d) %>% select(subgroup, est)
  ) %>% mutate(scenario = s, true_hr = exp(est)) %>% rename(true_loghr = est)
}) %>% bind_rows()

print(tidyr::pivot_wider(select(truth, subgroup, scenario, true_hr),
                         names_from = scenario, values_from = true_hr), n = Inf)
saveRDS(truth, "trial/outputs/truth.rds")
