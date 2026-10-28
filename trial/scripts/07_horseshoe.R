# Adds the horseshoe prior (M5) to the main study, on a subset: the first
# 200 trials of S1 (no subgroup differences) and S3 (a real ECOG difference)
# with 300 and 800 patients. These are exactly the same trials (same seeds)
# as in 03_study.R, so the horseshoe is compared with the standard forest
# plot, ridge and lasso trial by trial. Each fit is a Bayesian Cox model in
# Stan, which is why it isn't run on all 4,000 trials.
#
# Before running:
#   - scripts 02 and 03 (truth and raw results in trial/outputs/)
#   - bonsaiforest, brms and cmdstanr with CmdStan (cmdstanr::install_cmdstan();
#     Rtools on Windows)
#
# From the repo root:
#   Rscript trial/scripts/07_horseshoe.R
#
# Takes about 4 hours with 10 workers. A finished cell is skipped, so a
# stopped run can just be started again.
# Output: raw fits' estimates in trial/outputs/raw_hs/ (not in git), and
# horseshoe.csv, horseshoe_ecog.csv and fig_horseshoe.png in trial/results/.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(parallel)
})

n_rep <- 200
cells <- expand.grid(scenario = c("S1", "S3"), n = c(300, 800), stringsAsFactors = FALSE)
dir.create("trial/outputs/raw_hs", recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(10)
invisible(clusterCall(cl, function(p) .libPaths(p), .libPaths()))
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(dplyr); library(survival) })
  options(brms.backend = "cmdstanr")
  source("trial/R/simulate.R")
  source("trial/R/methods.R")
  base_fit <- NULL   # each worker compiles the model once, on its first trial
  one_rep <- function(scenario, n, rep) {
    seed <- 1e6 * match(scenario, c("S1", "S2", "S3", "S4")) + 1e3 * (n == 800) + rep   # as in 03
    d <- analysis_data(sim_trial(n, scenario, seed = seed))
    r <- tryCatch(m_horseshoe(d, base = base_fit, seed = rep), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    if (is.null(base_fit)) base_fit <<- r$fit
    list(est = r$est, max_rhat = r$max_rhat, divergent = r$divergent)
  }
}))

for (i in seq_len(nrow(cells))) {
  s <- cells$scenario[i]; n <- cells$n[i]
  f <- sprintf("trial/outputs/raw_hs/%s_n%d.rds", s, n)
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- parLapply(cl, seq_len(n_rep), function(r, s, n) one_rep(s, n, r), s = s, n = n)
  saveRDS(res, f)
  cat(s, n, "done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
stopCluster(cl)

# ---- compare with the other methods on the same trials ----
truth <- readRDS("trial/outputs/truth.rds") %>% select(scenario, subgroup, true_loghr)
mcse_mean <- function(x) sd(x) / sqrt(length(x))

est <- list(); diag <- list()
for (i in seq_len(nrow(cells))) {
  s <- cells$scenario[i]; n <- cells$n[i]
  hs <- readRDS(sprintf("trial/outputs/raw_hs/%s_n%d.rds", s, n))
  main <- readRDS(sprintf("trial/outputs/raw/%s_n%d.rds", s, n))
  for (r in seq_len(n_rep)) {
    if (is.null(hs[[r]])) next   # failed fit: that trial is left out for every method
    x <- main[[r]]
    est[[length(est) + 1]] <- bind_rows(
      mutate(x$naive, method = "Naive"),
      if (!is.null(x$ridge)) mutate(x$ridge, method = "Ridge"),
      if (!is.null(x$lasso)) mutate(x$lasso, method = "Lasso"),
      mutate(hs[[r]]$est, method = "Horseshoe")
    ) %>% mutate(scenario = s, n = n, rep = r)
    diag[[length(diag) + 1]] <- tibble(scenario = s, n = n, rep = r,
                                       max_rhat = hs[[r]]$max_rhat, divergent = hs[[r]]$divergent)
  }
}
est <- bind_rows(est) %>% left_join(truth, by = c("scenario", "subgroup"))
diag <- bind_rows(diag)

# RMSE as in 04_summarise.R: per-trial MSE first, MCSE by the delta method
t_rmse <- est %>%
  group_by(scenario, n, method, rep) %>%
  summarise(mse = mean((est - true_loghr)^2), .groups = "drop") %>%
  group_by(scenario, n, method) %>%
  summarise(trials = n(), rmse = sqrt(mean(mse)),
            rmse_mcse = mcse_mean(mse) / (2 * sqrt(mean(mse))), .groups = "drop")
t_diag <- diag %>% group_by(scenario, n) %>%
  summarise(fits = n(), share_rhat_over_1.01 = mean(max_rhat > 1.01),
            share_with_divergences = mean(divergent > 0), .groups = "drop")
res <- left_join(t_rmse, t_diag, by = c("scenario", "n")) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
print(as.data.frame(res))
write.csv(res, "trial/results/horseshoe.csv", row.names = FALSE)

# S3: does the method keep the real ECOG difference?
t_ecog <- est %>%
  filter(scenario == "S3", subgroup %in% c("ecog0", "ecog1")) %>%
  group_by(n, method, subgroup) %>%
  summarise(mean_hr = round(exp(mean(est)), 3), true_hr = round(exp(first(true_loghr)), 3),
            .groups = "drop")
print(t_ecog, n = Inf)
write.csv(t_ecog, "trial/results/horseshoe_ecog.csv", row.names = FALSE)

p <- res %>%
  mutate(method = factor(method, c("Naive", "Ridge", "Lasso", "Horseshoe")),
         scenario = recode(scenario, S1 = "S1 No effect", S3 = "S3 Effect only in ECOG 0"),
         size = paste("N =", n)) %>%
  ggplot(aes(method, rmse, fill = method)) + geom_col() +
  geom_errorbar(aes(ymin = rmse - 1.96 * rmse_mcse, ymax = rmse + 1.96 * rmse_mcse), width = 0.2) +
  facet_grid(size ~ scenario) +
  labs(x = NULL, y = "RMSE of subgroup log HR (lower is better)",
       title = "Adding the horseshoe prior",
       subtitle = "The same 200 trials per panel for every method; error bars = 95% Monte Carlo interval") +
  theme_bw() + theme(legend.position = "none")
ggsave("trial/results/fig_horseshoe.png", p, width = 8, height = 5.5, dpi = 150)
