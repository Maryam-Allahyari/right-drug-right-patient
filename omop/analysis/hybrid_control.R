# Can patients from outside the trial be added to its control arm?
# A hybrid control arm uses external patients to supplement the trial's own
# controls. Bayesian dynamic borrowing (psborrow2) decides how much to
# borrow from how similar the external controls are to the trial's controls:
# the commensurate prior borrows a lot when they agree and little when they
# don't. Full borrowing treats them as if they were trial patients.
#
# The trial is the showcase trial (S4, 800 patients, true overall HR about 0.72).
# Two sets of external controls:
#   - the routine-care cohort from OMOP (stage I, lives much longer: should
#     NOT be borrowed)
#   - the control arm of a second simulated trial with the same design
#     (exchangeable with the trial's controls: SHOULD be borrowed)
# OS is modelled as exponential, adjusted for age group and sex, the two
# covariates both sources have.
#
# Before running:
#   - trial/scripts/05_showcase.R and omop/analysis/compare_os.R run
#   - psborrow2, cmdstanr and CmdStan (cmdstanr::install_cmdstan(); on
#     Windows this needs Rtools), posterior, dplyr, ggplot2 installed
#
# From the repo root:
#   Rscript omop/analysis/hybrid_control.R
#
# Output in omop/results/: hybrid_control.csv and f_hybrid_control.png.
# Takes a few minutes, most of it compiling the Stan models.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(psborrow2)
})
source("trial/R/simulate.R")   # sim_trial(), analysis_data(), month
# Stan models are compiled by CmdStan through cmdstanr; stop early if it's missing
if (!requireNamespace("cmdstanr", quietly = TRUE) || is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  stop("needs cmdstanr and CmdStan: install.packages('cmdstanr', repos = 'https://stan-dev.r-universe.dev'); cmdstanr::install_cmdstan()")

out <- "omop/results"
trial_follow_up <- 30   # months; external follow-up is cut here too

# ---- the trial: both arms ----
trial <- read.csv("trial/results/showcase/adtte_os.csv") %>%
  transmute(source = "trial", trt = as.integer(ARM == "Drug X"), ext = 0L,
            time = AVALM, cens = CNSR,
            age65 = as.integer(AGEGR1 == ">=65"), male = as.integer(SEX == "M"))

# ---- external controls 1: routine care (OMOP cohort) ----
rwd <- read.csv("omop/results/adtte_os_rwd.csv") %>%
  transmute(source = "routine care", trt = 0L, ext = 1L,
            time = AVAL / month, cens = CNSR,
            age65 = as.integer(AGEGR1 == ">=65"), male = as.integer(SEX == "M"))

# ---- external controls 2: control arm of a second trial, same design ----
second <- analysis_data(sim_trial(800, "S4", seed = 20261028)) %>%
  filter(arm == 0) %>%
  transmute(source = "second trial", trt = 0L, ext = 1L,
            time = time, cens = 1L - event,
            age65 = as.integer(agegr == ">=65"), male = as.integer(sex == "M"))

# same follow-up window for everyone: alive at 30 months = censored there
cut_fu <- function(d) mutate(d, cens = if_else(time > trial_follow_up, 1L, cens),
                             time = pmin(time, trial_follow_up))

# ---- one analysis = data + outcome model + borrowing method ----
# Vague priors: normal(0, 1000) on the log baseline hazard and the log HRs.
# The commensurate prior's precision tau gets the usual vague gamma(0.001, 0.001):
# large tau = the two control groups look alike = borrow more.
fit_one <- function(ext_data, borrowing, label, seed) {
  dm <- create_data_matrix(cut_fu(bind_rows(trial, ext_data)),
                           outcome = c("time", "cens"), trt_flag_col = "trt",
                           ext_flag_col = "ext", covariates = c("age65", "male"))
  anls <- create_analysis_obj(
    data_matrix = dm,
    outcome = outcome_surv_exponential("time", "cens", baseline_prior = prior_normal(0, 1000)),
    borrowing = borrowing,
    treatment = treatment_details("trt", prior_normal(0, 1000)),
    covariates = add_covariates(c("age65", "male"), prior_normal(0, 1000)),
    quiet = TRUE
  )
  fit <- mcmc_sample(anls, iter_warmup = 1000, iter_sampling = 2000, chains = 4,
                     seed = seed)
  draws <- posterior::as_draws_df(fit$draws())
  hr_draws <- exp(draws$beta_trt)   # posterior HR, one value per draw
  tibble(analysis = label,
         hr = median(hr_draws), lower = quantile(hr_draws, .025),
         upper = quantile(hr_draws, .975),
         prob_benefit = mean(hr_draws < 1),
         tau = if ("tau" %in% names(draws)) median(draws$tau) else NA_real_,
         max_rhat = max(fit$summary()$rhat, na.rm = TRUE))
}

res <- bind_rows(
  fit_one(rwd, borrowing_none("ext"), "Trial only (no borrowing)", 1),
  fit_one(rwd, borrowing_full("ext"), "Routine care, full borrowing", 2),
  fit_one(rwd, borrowing_hierarchical_commensurate("ext", prior_gamma(0.001, 0.001)),
          "Routine care, dynamic borrowing", 3),
  fit_one(second, borrowing_full("ext"), "Second trial, full borrowing", 4),
  fit_one(second, borrowing_hierarchical_commensurate("ext", prior_gamma(0.001, 0.001)),
          "Second trial, dynamic borrowing", 5)
) %>%
  mutate(ci_width = upper - lower, across(where(is.numeric), ~ signif(.x, 3)))
print(res)
write.csv(res, file.path(out, "hybrid_control.csv"), row.names = FALSE)

# ---- figure: HR with 95% credible interval per analysis ----
p <- ggplot(res, aes(hr, factor(analysis, levels = rev(analysis)))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2, 3)) +
  labs(x = "HR for OS, Drug X vs chemotherapy (log scale), with 95% credible interval",
       y = NULL, title = "Adding external controls to the trial's control arm",
       subtitle = "Showcase trial (800 patients); external controls cut at 30 months") +
  theme_bw()
ggsave(file.path(out, "f_hybrid_control.png"), p, width = 8, height = 4, dpi = 150)
