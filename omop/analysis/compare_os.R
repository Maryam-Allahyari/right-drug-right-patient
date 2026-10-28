# Overall survival (OS) of the routine-care cohort next to the trial's control
# arm. Both get the same treatment (cisplatin + paclitaxel) and the same
# endpoint, built the same way: the cohort is turned into ADSL and ADTTE
# with admiral, exactly like the trial, and then analysed with the same
# tern tables.
#
# Before running:
#   - the cohort built with omop/cohort/build_cohort.R
#   - the trial showcase run (trial/scripts/05_showcase.R), for adtte_os.csv
#   - DatabaseConnector, SqlRender, config, admiral, survival, ggsurvfit,
#     ggplot2, rtables, tern and dplyr installed
#
# From the repo root:
#   Rscript omop/analysis/compare_os.R
#
# Output in omop/results/: the cohort's ADSL and ADTTE as csv, the baseline
# and OS tables as text, the hazard ratios as csv and the Kaplan-Meier plot.

suppressPackageStartupMessages({
  library(dplyr); library(admiral); library(survival)
  library(ggplot2); library(ggsurvfit); library(rtables); library(tern)
})

out <- "omop/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
month <- 30.4375   # same as the trial
trial_follow_up <- 30   # months, the longest follow-up in the simulated trial

# ---- the cohort from the database ----
cfg <- config::get(file = "omop/config.yml")
cd <- DatabaseConnector::createConnectionDetails(
  dbms = cfg$dbms, server = cfg$server, port = cfg$port,
  user = cfg$user, password = cfg$password,
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)
con <- DatabaseConnector::connect(cd)
# Stage comes from the source code: Synthea codes the diagnosis as e.g.
# "Non-small cell carcinoma of lung, TNM stage 1", but that maps to the plain
# "Non-small cell lung cancer" concept, so the stage is only left in
# condition_source_value. One stage code per patient.
sql <- SqlRender::render("
  SELECT c.subject_id, c.cohort_start_date, c.cohort_end_date, cand.age,
         p.gender_source_value AS sex, d.death_date, s.stage
  FROM @res.cohort c
  JOIN @res.nsclc_candidates cand ON cand.person_id = c.subject_id
  JOIN @cdm.person p ON p.person_id = c.subject_id
  LEFT JOIN @cdm.death d ON d.person_id = c.subject_id
  LEFT JOIN (
    SELECT person_id, MAX(CASE condition_source_value
             WHEN '424132000' THEN 'I'  WHEN '425048006' THEN 'II'
             WHEN '422968005' THEN 'III' WHEN '423121009' THEN 'IV' END) AS stage
    FROM @cdm.condition_occurrence
    WHERE condition_source_value IN ('424132000', '425048006', '422968005', '423121009')
    GROUP BY person_id) s ON s.person_id = c.subject_id
  WHERE c.cohort_definition_id = 1", cdm = cfg$cdm_schema, res = cfg$results_schema)
coh <- as_tibble(DatabaseConnector::querySql(con, SqlRender::translate(sql, cfg$dbms)))
DatabaseConnector::disconnect(con)
names(coh) <- toupper(names(coh))

# ---- ADSL ----
# The first platinum dose takes the place of randomisation (TRTSDT is the
# start date for OS). Patients are last known alive at the end of their
# records; a death after that couldn't be seen, the cohort check found none.
adsl <- coh %>%
  transmute(
    STUDYID = "RWD-NSCLC",
    USUBJID = paste0("RWD-", SUBJECT_ID),
    ARM = "Platinum chemotherapy (routine care)",
    AGE = as.numeric(AGE), AGEU = "YEARS",
    AGEGR1 = if_else(AGE >= 65, ">=65", "<65"),
    SEX,
    STAGE = coalesce(STAGE, "Unknown"),
    TRTSDT = as.Date(COHORT_START_DATE),
    DTHDT = as.Date(DEATH_DATE),
    DTHFL = if_else(!is.na(DTHDT), "Y", NA_character_),
    LSTALVDT = coalesce(DTHDT, as.Date(COHORT_END_DATE))
  )
stopifnot(all(is.na(adsl$DTHDT) | adsl$DTHDT >= adsl$TRTSDT))

# ---- ADTTE: OS with the same admiral event and censoring rules as the trial ----
adtte <- derive_param_tte(
  dataset_adsl = adsl, start_date = TRTSDT,
  event_conditions = list(death_event),
  censor_conditions = list(lastalive_censor),
  source_datasets = list(adsl = adsl),
  set_values_to = exprs(PARAMCD = "OS", PARAM = "Overall Survival")
) %>%
  # AVAL in days, counting the start day as day 1, like the trial
  derive_vars_duration(new_var = AVAL, start_date = STARTDT, end_date = ADT) %>%
  derive_vars_merged(dataset_add = adsl, by_vars = exprs(STUDYID, USUBJID),
                     new_vars = exprs(ARM, AGE, AGEGR1, SEX, STAGE))

write.csv(adsl, file.path(out, "adsl_rwd.csv"), row.names = FALSE)
write.csv(adtte, file.path(out, "adtte_os_rwd.csv"), row.names = FALSE)

# ---- both groups in one dataset ----
sources <- c("Trial control arm", "Routine care")   # trial first, it's the reference
trial_ctrl <- read.csv("trial/results/showcase/adtte_os.csv") %>%
  filter(ARM == "Platinum chemotherapy")
# the trial has no stage variable; its patients have stage IV disease by design (SAP)
os <- bind_rows(
  transmute(trial_ctrl, USUBJID, SOURCE = sources[1], AGE, AGEGR1, SEX,
            STAGE = "IV (trial eligibility)", AVAL, CNSR),
  transmute(adtte, USUBJID, SOURCE = sources[2], AGE, AGEGR1, SEX, STAGE, AVAL, CNSR)
) %>%
  mutate(SOURCE = factor(SOURCE, levels = sources),
         AGEGR1 = factor(AGEGR1, levels = c("<65", ">=65")), SEX = factor(SEX),
         STAGE = factor(STAGE), AVALM = AVAL / month, is_event = CNSR == 0)

# ---- baseline: age and sex are all the two sources have in common, plus stage ----
t_base <- basic_table(show_colcounts = TRUE,
                      title = "Baseline characteristics: trial control arm and routine care") %>%
  split_cols_by("SOURCE") %>%
  analyze_vars(c("AGE", "AGEGR1", "SEX", "STAGE")) %>%
  build_table(os)
writeLines(toString(t_base), file.path(out, "t_baseline.txt"))
print(t_base)

# ---- OS table: events, KM median, OS at 12 and 24 months, unadjusted HR ----
t_os <- basic_table(show_colcounts = TRUE,
                    title = "Overall survival: routine care compared with the trial control arm") %>%
  split_cols_by("SOURCE", ref_group = sources[1]) %>%
  count_values("is_event", values = TRUE, .labels = c(count_fraction = "Patients who died")) %>%
  surv_time(vars = "AVALM", var_labels = "OS (months)", is_event = "is_event",
            table_names = "time") %>%
  surv_timepoint(vars = "AVALM", var_labels = "Months", is_event = "is_event",
                 time_point = c(12, 24), method = "surv", table_names_suffix = "_tp",
                 .labels = c(event_free_rate = "Alive (%)")) %>%
  coxph_pairwise(vars = "AVALM", is_event = "is_event", var_labels = "Unadjusted analysis",
                 table_names = "cox") %>%
  build_table(os)
writeLines(toString(t_os), file.path(out, "t_os.txt"))
print(t_os)

# ---- hazard ratios over the trial's follow-up window ----
# The cohort is followed for up to 6 years, the trial for 30 months. For the
# HR both are cut at 30 months (alive patients censored there), so the two
# are compared over the same period. Adjusting for age group and sex shows
# whether those two differences explain the gap; stage, ECOG and histology
# aren't in the routine-care data, so they can't be adjusted for.
os30 <- os %>%
  mutate(is_event = is_event & AVALM <= trial_follow_up,
         AVALM = pmin(AVALM, trial_follow_up))
hr <- function(f, label) {
  fit <- coxph(f, data = os30)
  ci <- exp(confint(fit)["SOURCERoutine care", ])
  tibble(model = label, hr = exp(coef(fit)[["SOURCERoutine care"]]),
         lower = ci[[1]], upper = ci[[2]])
}
hrs <- bind_rows(
  hr(Surv(AVALM, is_event) ~ SOURCE, "Unadjusted, first 30 months"),
  hr(Surv(AVALM, is_event) ~ SOURCE + AGEGR1 + SEX, "Adjusted for age group and sex, first 30 months")
) %>% mutate(across(hr:upper, ~ round(.x, 2)))
print(hrs)
write.csv(hrs, file.path(out, "hazard_ratios.csv"), row.names = FALSE)

# ---- Kaplan-Meier plot ----
km <- survfit2(Surv(AVALM, is_event) ~ SOURCE, data = os) %>%
  ggsurvfit() + add_confidence_interval() + add_risktable() +
  geom_vline(xintercept = trial_follow_up, linetype = 3, colour = "grey40") +
  scale_x_continuous(breaks = seq(0, 72, 12)) +
  labs(x = "Months from randomisation (trial) or first platinum dose (routine care)", y = "Overall survival",
       title = "Overall survival: trial control arm and routine care",
       subtitle = "Both treated with cisplatin + paclitaxel. Dotted line: end of trial follow-up")
ggsave(file.path(out, "f_km_os.png"), km, width = 8, height = 6, dpi = 150)
