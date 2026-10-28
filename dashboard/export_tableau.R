# Exports the project's results as tidy CSV files for a Tableau dashboard
# (dashboard/README.md says how the dashboards are built). Long format, one
# row per number, so Tableau can filter and colour by any column. Nothing
# is re-estimated here except the Kaplan-Meier curves and the showcase
# forest plot, which are recomputed from the ADaM data so Tableau can draw them.
#
# Before running: the scripts in trial/scripts and omop/ (their results in
# trial/results and omop/results).
#
# From the repo root:
#   Rscript dashboard/export_tableau.R
#
# Output: dashboard/data/*.csv

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival) })

out <- "dashboard/data"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
month <- 30.4375
save <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "")
scen_lab <- c(S1 = "S1 No effect", S2 = "S2 Same effect for all",
              S3 = "S3 Effect only in ECOG 0", S4 = "S4 Continuous biomarker",
              S5 = "S5 Delayed effect")

# ---- Part 1: the simulation study, one row per measure ----
s <- read.csv("trial/results/summary.csv")
props <- c(power_overall = "Overall effect significant",
           any_interaction = "Any interaction test significant",
           ecog_power = "ECOG interaction found",
           other_vars_false = "Another variable wrongly flagged",
           any_misleading = "At least one misleading subgroup",
           bm_cutoff = "Biomarker effect found: cut-off",
           bm_linear = "Biomarker effect found: linear",
           bm_spline = "Biomarker effect found: spline",
           naive_coverage = "Coverage of 95% CIs")
measures <- bind_rows(
  lapply(names(props), function(m) s %>%
    transmute(scenario, n, measure = props[[m]], method = NA_character_,
              value = .data[[m]], mcse = .data[[paste0(m, "_mcse")]])),
  lapply(c("Naive", "Ridge", "Lasso"), function(m) s %>%
    transmute(scenario, n, measure = "RMSE of subgroup log HR", method = m,
              value = .data[[paste0("rmse_", m)]], mcse = .data[[paste0("rmse_mcse_", m)]]))
) %>% mutate(trials = 500)
# the horseshoe was run on a subset (200 trials of S1 and S3); kept apart by `trials`
if (file.exists("trial/results/horseshoe.csv")) {
  measures <- bind_rows(measures, read.csv("trial/results/horseshoe.csv") %>%
    transmute(scenario, n, measure = "RMSE of subgroup log HR", method,
              value = rmse, mcse = rmse_mcse, trials))
}
measures <- mutate(measures, scenario_label = scen_lab[scenario], .after = scenario)
save(measures, "study_measures.csv")

ecog <- read.csv("trial/results/s3_ecog_signal.csv") %>% mutate(trials = 500)
if (file.exists("trial/results/horseshoe_ecog.csv"))
  ecog <- bind_rows(ecog, read.csv("trial/results/horseshoe_ecog.csv") %>% mutate(trials = 200))
save(ecog, "s3_ecog_signal.csv")

save(read.csv("trial/results/delayed_effect.csv") %>%
       mutate(scenario_label = scen_lab[scenario], .after = scenario),
     "delayed_effect.csv")

# ---- Kaplan-Meier curves, as steps Tableau can draw ----
km_steps <- function(time, event, group, panel) {
  bind_rows(lapply(split(seq_along(time), group), function(i) {
    f <- survfit(Surv(time[i], event[i]) ~ 1)
    tibble(panel = panel, group = group[i][1],
           months = c(0, f$time), survival = c(1, f$surv),
           lower = c(1, f$lower), upper = c(1, f$upper),
           at_risk = c(length(i), f$n.risk))
  }))
}
adtte <- readRDS("trial/results/showcase/adtte.rds")
os <- filter(adtte, PARAMCD == "OS")
trial_ctrl <- read.csv("trial/results/showcase/adtte_os.csv") %>% filter(ARM == "Platinum chemotherapy")
rwd <- read.csv("omop/results/adtte_os_rwd.csv")
km <- bind_rows(
  km_steps(os$AVAL / month, os$CNSR == 0, as.character(os$ARM), "Showcase trial: OS by arm"),
  km_steps(c(trial_ctrl$AVALM, rwd$AVAL / month), c(trial_ctrl$CNSR, rwd$CNSR) == 0,
           rep(c("Trial control arm", "Routine care"), c(nrow(trial_ctrl), nrow(rwd))),
           "Trial control arm vs routine care")
)
save(km, "km_curves.csv")

# ---- forest plot of the showcase trial (OS), one row per subgroup ----
os <- mutate(os, ARM = relevel(factor(as.character(ARM)), ref = "Platinum chemotherapy"))
sg <- c(Sex = "SEX", `Age group` = "AGEGR1", Region = "REGION1", ECOG = "ECOG",
        Histology = "HISTOL", `Biomarker group` = "BMRKGR1")
forest <- bind_rows(
  { fit <- coxph(Surv(AVAL, CNSR == 0) ~ ARM, data = os)
    tibble(variable = "All patients", level = "All patients", patients = nrow(os),
           events = sum(os$CNSR == 0), hr = exp(coef(fit)), lower = exp(confint(fit))[1],
           upper = exp(confint(fit))[2]) },
  lapply(names(sg), function(v) bind_rows(lapply(sort(unique(as.character(os[[sg[[v]]]]))), function(l) {
    dd <- os[as.character(os[[sg[[v]]]]) == l, ]
    fit <- coxph(Surv(AVAL, CNSR == 0) ~ ARM, data = dd)
    tibble(variable = v, level = l, patients = nrow(dd), events = sum(dd$CNSR == 0),
           hr = exp(coef(fit)), lower = exp(confint(fit))[1], upper = exp(confint(fit))[2])
  })))
) %>% mutate(across(hr:upper, ~ round(.x, 3)))
save(forest, "forest_showcase.csv")

# ---- Part 2: the real-world cohort ----
baseline <- bind_rows(
  tibble(source = "Trial control arm", age = trial_ctrl$AGE, sex = trial_ctrl$SEX),
  tibble(source = "Routine care", age = rwd$AGE, sex = rwd$SEX)
) %>%
  group_by(source) %>%
  summarise(patients = n(), `Median age` = median(age), `Aged 65 or over (%)` = 100 * mean(age >= 65),
            `Male (%)` = 100 * mean(sex == "M"), .groups = "drop") %>%
  pivot_longer(-source, names_to = "measure", values_to = "value") %>%
  mutate(value = round(value, 1))
save(baseline, "baseline_trial_vs_routine_care.csv")
save(read.csv("omop/cohort/attrition.csv"), "cohort_attrition.csv")
save(read.csv("omop/results/hybrid_control.csv"), "hybrid_control.csv")
save(read.csv("omop/quality/dqd_overview.csv", check.names = FALSE), "data_quality_overview.csv")

cat("written:", paste(list.files(out), collapse = ", "), "\n")
