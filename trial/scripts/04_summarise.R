# Turns the raw study results into the performance measures from the plan,
# each with its Monte Carlo standard error (MCSE), plus three figures.
#
# Before running: scripts 02 and 03 (they write trial/outputs/).
#
# From the repo root:
#   Rscript trial/scripts/04_summarise.R
#
# Output in trial/results/: summary.csv (one row per scenario and size),
# s3_ecog_signal.csv, fig_rmse.png, fig_s3_signal.png, fig_biomarker_power.png

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

raw <- "trial/outputs"
out <- "trial/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sg_vars <- c("sex", "agegr", "region", "ecog", "hist", "bmhigh")

truth <- readRDS(file.path(raw, "truth.rds")) %>% select(scenario, subgroup, true_loghr)

# MCSE of a proportion, and of the mean of per-trial values
mcse_p    <- function(p, n) sqrt(p * (1 - p) / n)
mcse_mean <- function(x) sd(x) / sqrt(length(x))

# ---- read and tidy the raw results ----
est <- list(); tests <- list(); overall <- list()
for (f in list.files(file.path(raw, "raw"), full.names = TRUE)) {
  cell <- regmatches(basename(f), regexec("(S[1-4])_n(\\d+)", basename(f)))[[1]]
  res <- readRDS(f)
  for (r in seq_along(res)) {
    x <- res[[r]]
    key <- tibble(scenario = cell[2], n = as.integer(cell[3]), rep = r)
    est[[length(est) + 1]] <- bind_rows(
      mutate(x$naive, method = "Naive"),
      if (!is.null(x$ridge)) mutate(x$ridge, method = "Ridge"),
      if (!is.null(x$lasso)) mutate(x$lasso, method = "Lasso")
    ) %>% bind_cols(key[rep(1, nrow(.)), ])
    tests[[length(tests) + 1]] <- bind_cols(
      key, as_tibble(as.list(x$interaction)),
      as_tibble(as.list(setNames(x$biomarker, paste0("bm_", names(x$biomarker))))))
    overall[[length(overall) + 1]] <- bind_cols(
      key, tibble(est = x$overall[["est"]], p = x$overall[["p"]], events = x$events))
  }
}
est <- bind_rows(est) %>% left_join(truth, by = c("scenario", "subgroup"))
tests <- bind_rows(tests); overall <- bind_rows(overall)

# any shrinkage fits that failed show up here as fewer than 500 trials
fits <- est %>% distinct(scenario, n, rep, method) %>% count(scenario, n, method, name = "fits")
print(filter(fits, fits < 500))

# ---- 1. the overall result of each trial ----
t_overall <- overall %>% group_by(scenario, n) %>%
  summarise(trials = n(), mean_events = mean(events), mean_hr = exp(mean(est)),
            power_overall = mean(p < 0.05), .groups = "drop")

# ---- 2. interaction tests: false alarms and power ----
t_het <- tests %>% group_by(scenario, n) %>%
  summarise(
    any_interaction = mean(apply(across(all_of(sg_vars)), 1, min) < 0.05),
    # in S3 only ECOG is real, so a hit on any of the other four
    # (not ECOG, not the biomarker) is a false alarm
    ecog_power = mean(ecog < 0.05),
    other_vars_false = mean(apply(across(all_of(setdiff(sg_vars, c("ecog", "bmhigh")))), 1, min) < 0.05),
    bm_cutoff = mean(bm_cutoff < 0.05), bm_linear = mean(bm_linear < 0.05),
    bm_spline = mean(bm_spline < 0.05),
    .groups = "drop")

# ---- 3. misleading subgroups in the naive forest plot ----
# "misleading" = CI excludes 1 on the wrong side of the truth, or excludes 1
# when the true HR is (close to) 1
t_mislead <- est %>% filter(method == "Naive") %>%
  mutate(wrong = (lo > 0 & true_loghr < -0.05) | (hi < 0 & true_loghr > 0.05) |
                 (abs(true_loghr) <= 0.05 & (hi < 0 | lo > 0))) %>%
  group_by(scenario, n, rep) %>% summarise(any_wrong = any(wrong), .groups = "drop") %>%
  group_by(scenario, n) %>% summarise(any_misleading = mean(any_wrong), .groups = "drop")

# ---- 4. accuracy (RMSE) and coverage ----
# per-trial mean squared error first, so the MCSE comes from the spread
# between trials; MCSE of RMSE by the delta method: MCSE(MSE) / (2 * RMSE)
t_rmse <- est %>% group_by(scenario, n, method, rep) %>%
  summarise(mse = mean((est - true_loghr)^2), .groups = "drop") %>%
  group_by(scenario, n, method) %>%
  summarise(rmse = sqrt(mean(mse)), rmse_mcse = mcse_mean(mse) / (2 * sqrt(mean(mse))),
            .groups = "drop")
t_cov <- est %>% filter(method == "Naive") %>% group_by(scenario, n, rep) %>%
  summarise(cov = mean(lo <= true_loghr & true_loghr <= hi), .groups = "drop") %>%
  group_by(scenario, n) %>%
  summarise(naive_coverage = mean(cov), naive_coverage_mcse = mcse_mean(cov), .groups = "drop")

# ---- 5. S3: does shrinkage keep the real ECOG difference? ----
t_ecog <- est %>% filter(scenario == "S3", subgroup %in% c("ecog0", "ecog1")) %>%
  group_by(n, method, subgroup) %>%
  summarise(mean_hr = exp(mean(est)), true_hr = exp(first(true_loghr)), .groups = "drop")

# one wide table, proportions with their MCSE next to them
props <- c("power_overall", "any_interaction", "ecog_power", "other_vars_false",
           "bm_cutoff", "bm_linear", "bm_spline", "any_misleading")
summary_tbl <- t_overall %>%
  left_join(t_het, by = c("scenario", "n")) %>%
  left_join(t_mislead, by = c("scenario", "n")) %>%
  left_join(t_cov, by = c("scenario", "n")) %>%
  left_join(pivot_wider(t_rmse, names_from = method, values_from = c(rmse, rmse_mcse)),
            by = c("scenario", "n")) %>%
  mutate(across(all_of(props), ~ mcse_p(.x, trials), .names = "{.col}_mcse"))

print(as.data.frame(summary_tbl), digits = 3)
print(t_ecog, n = Inf)
write.csv(summary_tbl, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(t_ecog, file.path(out, "s3_ecog_signal.csv"), row.names = FALSE)

# ---- figures ----
scen_lab <- c(S1 = "S1 No effect", S2 = "S2 Same effect for all",
              S3 = "S3 Effect only in ECOG 0", S4 = "S4 Continuous biomarker")
methods <- c("Naive", "Ridge", "Lasso")

p1 <- t_rmse %>%
  mutate(method = factor(method, methods), size = paste("N =", n)) %>%
  ggplot(aes(method, rmse, fill = method)) + geom_col() +
  geom_errorbar(aes(ymin = rmse - 1.96 * rmse_mcse, ymax = rmse + 1.96 * rmse_mcse), width = 0.2) +
  facet_grid(size ~ scenario, labeller = labeller(scenario = scen_lab)) +
  labs(x = NULL, y = "RMSE of subgroup log HR (lower is better)",
       title = "How far the subgroup estimates are from the truth (OS)",
       subtitle = "500 trials per bar; error bars = 95% Monte Carlo interval") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out, "fig_rmse.png"), p1, width = 10, height = 6, dpi = 150)

true_s3 <- truth %>% filter(scenario == "S3", subgroup %in% c("ecog0", "ecog1")) %>%
  mutate(subgroup = recode(subgroup, ecog0 = "ECOG 0", ecog1 = "ECOG 1"))
p2 <- est %>% filter(scenario == "S3", subgroup %in% c("ecog0", "ecog1")) %>%
  mutate(method = factor(method, methods), size = paste("N =", n),
         subgroup = recode(subgroup, ecog0 = "ECOG 0", ecog1 = "ECOG 1")) %>%
  ggplot(aes(method, exp(est), fill = subgroup)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_hline(data = true_s3, aes(yintercept = exp(true_loghr), colour = subgroup),
             linetype = 2, show.legend = FALSE) +
  scale_fill_manual(values = c("ECOG 0" = "#1b9e77", "ECOG 1" = "#d95f02")) +
  scale_colour_manual(values = c("ECOG 0" = "#1b9e77", "ECOG 1" = "#d95f02")) +
  facet_wrap(~ size) + scale_y_log10() +
  labs(x = NULL, y = "Estimated HR for OS (log scale)", fill = NULL,
       title = "S3: a real subgroup effect (dashed lines = truth)",
       subtitle = "Shrinkage pulls the two groups towards each other") +
  theme_bw()
ggsave(file.path(out, "fig_s3_signal.png"), p2, width = 9, height = 5, dpi = 150)

p3 <- t_het %>% select(scenario, n, bm_cutoff, bm_linear, bm_spline) %>%
  pivot_longer(starts_with("bm_"), names_to = "model", values_to = "share") %>%
  mutate(model = recode(model, bm_cutoff = "Cut-off (High/Low)", bm_linear = "Linear (log)",
                        bm_spline = "Spline (log)"), size = paste("N =", n)) %>%
  ggplot(aes(scenario, share, fill = model)) + geom_col(position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = 2) + facet_wrap(~ size) +
  scale_x_discrete(labels = scen_lab) +
  labs(x = NULL, y = "Share of trials with interaction p < 0.05", fill = NULL,
       title = "Biomarker x treatment interaction: false alarms (S1-S3) and power (S4)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
ggsave(file.path(out, "fig_biomarker_power.png"), p3, width = 10, height = 5, dpi = 150)
