# Simulates one first-line NSCLC trial: Drug X vs platinum chemotherapy.
# The design (covariates, scenarios, follow-up) is in docs/sap.md.
#
# sim_trial()     patient-level data, including the scan schedule
# pfs_os_truth()  OS and PFS exactly as the plan defines them; the ADaM
#                 pipeline has to reproduce these values
# analysis_data() one row per patient for the big study (OS, covariates)
#
# All times are in days from randomisation.

library(dplyr)

month <- 30.4375

# true log HR for each patient, by scenario
scenario_loghr <- function(scenario, d) {
  switch(scenario,
    S1 = rep(0, nrow(d)),
    S2 = rep(log(0.70), nrow(d)),
    S3 = ifelse(d$ecog == "0", log(0.55), 0),
    S4 = log(0.70) - 0.5 * log(d$bmrk / 5)
  )
}

# countries per region, so DM has a real COUNTRY and ADSL can derive the region
countries <- list("North America" = c("USA", "CAN"),
                  "Europe"        = c("DEU", "FRA", "ESP", "GBR"),
                  "Asia"          = c("CHN", "JPN", "KOR"))

sim_trial <- function(n, scenario, seed = NULL, max_scans = 24) {
  if (!is.null(seed)) set.seed(seed)
  d <- data.frame(
    id     = seq_len(n),
    arm    = sample(rep(0:1, length.out = n)),     # 1:1, exactly balanced
    sex    = sample(c("F", "M"), n, replace = TRUE, prob = c(.4, .6)),
    age    = pmin(pmax(round(rnorm(n, 63, 9)), 30), 85),
    region = sample(names(countries), n, replace = TRUE, prob = c(.3, .4, .3)),
    ecog   = sample(c("0", "1"), n, replace = TRUE, prob = c(.4, .6)),
    hist   = sample(c("Non-squamous", "Squamous"), n, replace = TRUE, prob = c(.7, .3)),
    bmrk   = round(exp(rnorm(n, log(5), 0.6)), 2)
  )
  d$country <- vapply(d$region, function(r) sample(countries[[r]], 1), "")
  d$agegr  <- ifelse(d$age >= 65, ">=65", "<65")
  d$bmhigh <- ifelse(d$bmrk > 5, "High", "Low")
  d$loghr  <- scenario_loghr(scenario, d)

  # Death and progression are two competing exponential times. The same
  # prognostic and treatment effects multiply both hazards, so OS and PFS
  # share one true HR. Medians are for the best-prognosis group; they give
  # about 11 months OS and 5 months PFS for the control arm as a whole.
  lp <- log(1.2) * (d$agegr == ">=65") + log(1.5) * (d$ecog == "1") +
    log(1.2) * (d$hist == "Squamous") + d$arm * d$loghr
  t_death <- rexp(n, log(2) / (16 * month) * exp(lp))
  t_prog  <- rexp(n, log(2) / (12.5 * month) * exp(lp))
  d$t_death <- ceiling(t_death)
  d$t_prog  <- ceiling(t_prog)     # only seen if it comes before death

  d$rand_day <- floor(runif(n, 0, 18 * month))     # day of randomisation in the study
  cutoff     <- floor(30 * month) - d$rand_day      # data cut-off, 30 months after the first patient
  dropout    <- ceiling(rexp(n, -log(0.95) / 365.25))   # 5% lost per year
  d$fu_end   <- pmin(cutoff, dropout)
  d$dropped  <- dropout < cutoff

  # scan days: every 42 days, +/- 3 days
  scans <- outer(rep(1, n), 42 * seq_len(max_scans)) +
    matrix(sample(-3:3, n * max_scans, replace = TRUE), n)
  ok <- scans <= d$fu_end & scans < d$t_death
  pd_scan <- ok & scans >= d$t_prog
  first_pd <- apply(pd_scan, 1, function(x) if (any(x)) which(x)[1] else NA)
  # scans continue up to and including the one that shows progression
  n_scans <- ifelse(is.na(first_pd), rowSums(ok), first_pd)
  d$pd_day <- ifelse(is.na(first_pd), NA, scans[cbind(seq_len(n), pmax(first_pd, 1))])
  d$n_scans <- n_scans
  d$last_scan <- ifelse(n_scans > 0, scans[cbind(seq_len(n), pmax(n_scans, 1))], NA)

  # partial response from the 2nd scan on. 40% of patients are responders on
  # chemotherapy, but some progress or die before a 2nd scan, which gives an
  # overall response rate of about 25%. More respond on Drug X when it works.
  d$responder <- runif(n) < plogis(qlogis(0.40) - 1.5 * d$arm * d$loghr)

  attr(d, "scans") <- scans
  d
}

# OS and PFS by the plan's rules (the answer key for the ADaM pipeline)
pfs_os_truth <- function(d) {
  death_obs <- d$t_death <= d$fu_end
  trt_end <- pmin(ifelse(is.na(d$pd_day), Inf, d$pd_day),
                  ifelse(death_obs, d$t_death, Inf), d$fu_end)
  d %>% mutate(
    pfs_event = as.integer(!is.na(pd_day) | death_obs),
    pfs_day   = case_when(!is.na(pd_day) ~ pd_day,
                          death_obs ~ t_death,
                          !is.na(last_scan) ~ last_scan,
                          TRUE ~ 0),
    os_event  = as.integer(death_obs),
    # survival is followed up after progression too (SS domain), so anyone
    # alive is censored at the cut-off or when they were lost, not at the last scan
    os_day    = ifelse(death_obs, t_death, fu_end),
    trt_end   = trt_end,
    # ADaM counts the start day as day 1, hence the +1
    PFS_AVAL  = pfs_day + 1,
    OS_AVAL   = os_day + 1
  )
}

# one row per patient with OS in months, used by the big study
# (no SDTM round trip; scripts/01 shows the pipeline gives the same values)
analysis_data <- function(d) {
  d <- pfs_os_truth(d)
  data.frame(
    arm    = factor(d$arm, levels = 0:1),
    time   = d$OS_AVAL / month,
    event  = d$os_event,
    sex    = factor(d$sex),
    agegr  = factor(d$agegr, levels = c("<65", ">=65")),
    region = factor(d$region, levels = names(countries)),
    ecog   = factor(d$ecog),
    hist   = factor(d$hist, levels = c("Non-squamous", "Squamous")),
    bmhigh = factor(d$bmhigh, levels = c("Low", "High")),
    bmrk   = d$bmrk
  )
}
