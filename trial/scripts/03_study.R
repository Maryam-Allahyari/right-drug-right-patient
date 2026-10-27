# The simulation study: 4 scenarios x 2 trial sizes x 500 trials.
# Each trial gets the overall Cox model, the naive forest plot, the
# interaction tests, ridge, lasso and the three biomarker models.
#
# Before running:
#   install.packages(c("bonsaiforest", "dplyr", "survival"))
#
# From the repo root:
#   Rscript trial/scripts/03_study.R
#
# Output: one .rds file of raw results per scenario and size in
# trial/outputs/raw (not in git, rebuilt by this script).
# A finished file is skipped, so a stopped run can just be started again.

suppressPackageStartupMessages(library(parallel))

n_rep <- 500
cells <- expand.grid(scenario = c("S1", "S2", "S3", "S4"), n = c(300, 800),
                     stringsAsFactors = FALSE)
dir.create("trial/outputs/raw", recursive = TRUE, showWarnings = FALSE)

one_rep <- function(scenario, n, rep) {
  # every trial has its own fixed seed, so any single trial can be rerun
  seed <- 1e6 * match(scenario, c("S1", "S2", "S3", "S4")) + 1e3 * (n == 800) + rep
  d <- analysis_data(sim_trial(n, scenario, seed = seed))
  list(
    overall = overall_cox(d),
    naive = m_naive(d),
    interaction = m_interaction(d),
    # shrinkage can fail on an odd trial; keep going and count the failures later
    ridge = tryCatch(m_elastic(d, 0), error = function(e) NULL),
    lasso = tryCatch(m_elastic(d, 1), error = function(e) NULL),
    biomarker = m_biomarker(d),
    events = sum(d$event)
  )
}

cl <- makeCluster(min(16, detectCores() - 4))
invisible(clusterCall(cl, function(p) .libPaths(p), .libPaths()))   # same packages as here
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(dplyr))
  source("trial/R/simulate.R")
  source("trial/R/methods.R")
}))
clusterExport(cl, "one_rep")

for (i in seq_len(nrow(cells))) {
  s <- cells$scenario[i]; n <- cells$n[i]
  f <- sprintf("trial/outputs/raw/%s_n%d.rds", s, n)
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- parLapply(cl, seq_len(n_rep), function(r, s, n) one_rep(s, n, r), s = s, n = n)
  saveRDS(res, f)
  cat(s, n, "done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
stopCluster(cl)
