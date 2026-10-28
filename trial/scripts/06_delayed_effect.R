# When the treatment effect starts late, is the HR still the right summary?
# Immunotherapies often show no difference for the first months and a
# clear one after. The HR assumes the effect is the same at all times
# (proportional hazards), so it averages the early and late parts. The
# restricted mean survival time (RMST) asks something simpler: how many
# months, on average, does a patient live within the first 24 months?
# The difference between the arms is in months, and it doesn't assume
# anything about how the effect changes over time.
#
# Compared on S2 (HR 0.70 from the start) and S5 (no effect for 4 months,
# then HR 0.60), 500 trials each with 300 and 800 patients:
#   - power of the log-rank test and of the RMST difference
#   - the average HR, and how often a test of proportional hazards notices
#     that the effect changes over time
#   - bias and coverage of the RMST difference, against the truth from one
#     trial of 200,000 patients
#
# From the repo root:
#   Rscript trial/scripts/06_delayed_effect.R
#
# Output in trial/results/: delayed_effect.csv and fig_delayed_effect.png.

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(ggplot2); library(parallel)
})
source("trial/R/simulate.R")

out <- "trial/results"
tau <- 24   # months; RMST up to here (minimum follow-up is 12 months, most have 24+)

# RMST per arm and their difference (Drug X - chemotherapy), with a z-test
rmst_diff <- function(d) {
  tab <- summary(survfit(Surv(time, event) ~ arm, data = d), rmean = tau)$table
  diff <- tab[2, "rmean"] - tab[1, "rmean"]
  se <- sqrt(tab[1, "se(rmean)"]^2 + tab[2, "se(rmean)"]^2)
  c(diff = diff, se = se, p = 2 * pnorm(-abs(diff / se)))
}

one_trial <- function(scenario, n, rep) {
  # seeds kept apart from the main study's
  seed <- 7e6 + 1e5 * match(scenario, c("S2", "S5")) + 1e3 * (n == 800) + rep
  d <- analysis_data(sim_trial(n, scenario, seed = seed))
  fit <- coxph(Surv(time, event) ~ arm, data = d)
  c(loghr = coef(fit)[[1]],
    p_logrank = summary(fit)$sctest[["pvalue"]],   # the score test is the log-rank test
    p_ph = cox.zph(fit)$table["arm", "p"],
    rmst_diff(d))
}

# ---- truth: one very large trial per scenario ----
truth <- lapply(c("S2", "S5"), function(s) {
  d <- analysis_data(sim_trial(2e5, s, seed = 777))
  tibble(scenario = s, true_rmst_diff = rmst_diff(d)[["diff"]],
         average_hr = exp(coef(coxph(Surv(time, event) ~ arm, data = d))[[1]]))
}) %>% bind_rows()
print(truth)

# ---- 500 trials per scenario and size, in parallel ----
cl <- makeCluster(min(16, detectCores() - 4))
invisible(clusterCall(cl, function(p) .libPaths(p), .libPaths()))
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(dplyr); library(survival) })
  source("trial/R/simulate.R")
}))
clusterExport(cl, c("one_trial", "rmst_diff", "tau"))
jobs <- expand.grid(rep = 1:500, n = c(300, 800), scenario = c("S2", "S5"),
                    stringsAsFactors = FALSE)
raw <- parLapply(cl, seq_len(nrow(jobs)), function(i, jobs)
  one_trial(jobs$scenario[i], jobs$n[i], jobs$rep[i]), jobs = jobs)
stopCluster(cl)
raw <- bind_cols(jobs, bind_rows(lapply(raw, function(x) as.list(x))))

# ---- performance, with Monte Carlo standard errors ----
mcse_p <- function(p, k) sqrt(p * (1 - p) / k)
res <- raw %>%
  left_join(truth, by = "scenario") %>%
  group_by(scenario, n) %>%
  summarise(
    trials = n(),
    power_logrank = mean(p_logrank < 0.05),
    power_rmst = mean(p < 0.05),
    mean_hr = exp(mean(loghr)),
    ph_rejected = mean(p_ph < 0.05),
    true_rmst_diff = first(true_rmst_diff),
    mean_rmst_diff = mean(diff),
    rmst_bias = mean(diff - true_rmst_diff),
    rmst_coverage = mean(abs(diff - true_rmst_diff) < 1.96 * se),
    .groups = "drop"
  ) %>%
  mutate(mcse_power_logrank = mcse_p(power_logrank, trials),
         mcse_power_rmst = mcse_p(power_rmst, trials),
         mcse_coverage = mcse_p(rmst_coverage, trials),
         across(where(is.double), ~ round(.x, 3)))
print(as.data.frame(res))
write.csv(res, file.path(out, "delayed_effect.csv"), row.names = FALSE)

# ---- figure: what the two effects look like, and the power of each test ----
km_df <- lapply(c("S2", "S5"), function(s) {
  d <- analysis_data(sim_trial(20000, s, seed = 778))
  f <- survfit(Surv(time, event) ~ arm, data = d)
  data.frame(scenario = s, time = f$time, surv = f$surv,
             arm = rep(c("Chemotherapy", "Drug X"), f$strata))
}) %>% bind_rows() %>%
  mutate(scenario = recode(scenario, S2 = "S2: HR 0.70 from the start",
                           S5 = "S5: no effect for 4 months, then HR 0.60"))
p_km <- ggplot(km_df, aes(time, surv, colour = arm)) +
  geom_step() + facet_wrap(~scenario) +
  geom_vline(xintercept = tau, linetype = 3) +
  labs(x = "Months", y = "Overall survival", colour = NULL,
       title = "Two ways a drug can work (20,000 simulated patients each)",
       subtitle = "Dotted line: 24 months, where the RMST is restricted") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out, "fig_delayed_effect.png"), p_km, width = 9, height = 4.5, dpi = 150)
