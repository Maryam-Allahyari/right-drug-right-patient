# The analysis methods from the plan, each applied to one trial's
# analysis data (from analysis_data(): OS time in months, event, subgroups).
# Estimates are log hazard ratios, Drug X vs chemotherapy.

library(dplyr)
library(survival)
library(splines)

sg_vars <- c("sex", "agegr", "region", "ecog", "hist", "bmhigh")

# subgroup labels without spaces, so they match bonsaiforest's output
sg_label <- function(v, l) gsub("[ -]", "", paste0(v, l))

# M1: separate Cox model in each subgroup (the usual forest plot)
m_naive <- function(d) {
  lapply(sg_vars, function(v) lapply(levels(d[[v]]), function(l) {
    fit <- coxph(Surv(time, event) ~ arm, data = d[d[[v]] == l, ])
    ci <- confint(fit)
    tibble(subgroup = sg_label(v, l), est = coef(fit)[[1]], lo = ci[1], hi = ci[2])
  })) %>% bind_rows()
}

# M2: likelihood ratio test of treatment x variable, one per variable
m_interaction <- function(d) {
  sapply(sg_vars, function(v) {
    f0 <- as.formula(paste("Surv(time, event) ~ arm +", v))
    f1 <- as.formula(paste("Surv(time, event) ~ arm *", v))
    anova(coxph(f0, d), coxph(f1, d))[2, "Pr(>|Chi|)"]
  })
}

# M3 (alpha = 0) and M4 (alpha = 1): ridge and lasso shrinkage from bonsaiforest.
# It returns HRs, so they're logged to compare with the other methods.
m_elastic <- function(d, alpha) {
  fit <- bonsaiforest::elastic_net("time", "arm", sg_vars, sg_vars, d, "survival",
                                   alpha = alpha, status = "event")
  summary(fit)$estimates %>%
    transmute(subgroup = gsub("[ -]", "", subgroup), est = log(trt.estimate))
}

# B1-B3: does the treatment effect change with the biomarker?
# Same test three ways: the assay cut-off, log biomarker as a line, and a spline.
m_biomarker <- function(d) {
  d$lb <- log(d$bmrk)
  lrt <- function(f0, f1) anova(coxph(f0, d), coxph(f1, d))[2, "Pr(>|Chi|)"]
  c(
    cutoff = lrt(Surv(time, event) ~ arm + bmhigh, Surv(time, event) ~ arm * bmhigh),
    linear = lrt(Surv(time, event) ~ arm + lb,     Surv(time, event) ~ arm * lb),
    spline = lrt(Surv(time, event) ~ arm + ns(lb, df = 3), Surv(time, event) ~ arm * ns(lb, df = 3))
  )
}

# overall treatment effect (the trial's main result)
overall_cox <- function(d) {
  fit <- coxph(Surv(time, event) ~ arm, data = d)
  c(est = coef(fit)[[1]], p = summary(fit)$coefficients[1, "Pr(>|z|)"])
}
