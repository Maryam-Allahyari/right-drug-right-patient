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

# M5: horseshoe prior, bonsaiforest's Bayesian shrinkage model fitted with
# brms/Stan. bonsaiforest::horseshoe() compiles a new Stan model on every
# call; here the first fit is compiled once (base = NULL) and later fits
# reuse it with new data (base = that fit), which is the same model and
# saves a compile per trial. The code otherwise follows horseshoe().
# brms builds formulas from the subgroup level names, so levels like "<65"
# or "North America" are renamed first and mapped back at the end.
m_horseshoe <- function(d, base = NULL, seed = 1) {
  levels(d$agegr)  <- c("lt65", "ge65")
  levels(d$region) <- gsub(" ", "", levels(d$region))
  levels(d$hist)   <- gsub("-", "", levels(d$hist))
  prep <- bonsaiforest:::preprocess("arm", sg_vars, sg_vars, d)
  design <- cbind(prep$design_main, prep$design_ia)
  y <- data.frame(time = d$time, event = d$event)
  # baseline hazard: M-splines with knots at the quartiles of the times, and
  # boundaries one step outside the data. The step is the smallest gap between
  # distinct times: times are whole days, so tied times would give a step of 0
  # and put data on the boundary, where the spline basis breaks down
  srt <- sort(unique(y$time))
  step <- min(diff(c(0, srt)))
  lim <- c(max(min(y$time) - step, 0), max(y$time) + step)
  knots <- quantile(y$time, c(.25, .5, .75), names = FALSE)
  form <- as.formula(paste0("time | cens(1 - event) + bhaz(Boundary.knots = ", deparse(lim),
                            ", knots = ", deparse(knots), ", intercept = FALSE) ~ a + b"))
  form_a <- as.formula(paste("a ~ 0 +", paste(colnames(prep$design_main), collapse = " + ")))
  form_b <- as.formula(paste("b ~ 0 +", paste(colnames(prep$design_ia), collapse = " + ")))
  bform <- brms::bf(form, nl = TRUE) + brms::lf(form_a) + brms::lf(form_b)
  args <- list(chains = 4, cores = 1, iter = 2000, warmup = 1000, seed = seed, refresh = 0,
               control = list(adapt_delta = 0.95), silent = 2)
  fit <- if (is.null(base)) {
    # main effects get a wide normal prior, treatment interactions the horseshoe
    do.call(brms::brm, c(list(bform, data = cbind(design, y), family = brms::brmsfamily("cox"),
      prior = brms::prior_string("normal(0, 5)", class = "b", nlpar = "a") +
              brms::prior_string("horseshoe(1)", class = "b", nlpar = "b")), args))
  } else {
    do.call(update, c(list(base, formula. = bform, newdata = cbind(design, y),
                           recompile = FALSE), args))
  }
  res <- structure(list(fit = fit, model = "horseshoe", resptype = "survival", data = d,
                        design_matrix = design, design_dummy = prep$design_dummy, y = y,
                        subgr_names = prep$subgr_names), class = c("bonsaiforest", "horseshoe"))
  np <- brms::nuts_params(fit)
  est <- summary(res)$estimates %>%
    transmute(subgroup = recode(subgroup, agegrlt65 = "agegr<65", agegrge65 = "agegr>=65"),
              est = log(trt.estimate), lo = log(trt.low), hi = log(trt.high))
  list(est = est, fit = fit,
       max_rhat = max(brms::rhat(fit), na.rm = TRUE),
       divergent = sum(np$Value[np$Parameter == "divergent__"]))
}
