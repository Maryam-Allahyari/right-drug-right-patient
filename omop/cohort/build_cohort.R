# Builds the real-world comparator cohort: NSCLC patients starting platinum
# chemotherapy (cisplatin + paclitaxel) for the first time. The definition is
# in nsclc_platinum.sql; this script runs it, counts how many patients each
# inclusion rule removes (attrition) and summarises who is left.
#
# Before running:
#   - the OMOP database loaded (omop/etl/load_omop.R)
#   - DatabaseConnector, SqlRender, dplyr and config installed
#
# From the repo root:
#   Rscript omop/cohort/build_cohort.R
#
# Output: the cohort in results.cohort (cohort_definition_id 1), every
# candidate with its rule flags in results.nsclc_candidates, and
# omop/cohort/attrition.csv and omop/cohort/cohort_summary.csv.

suppressPackageStartupMessages(library(dplyr))

cohort_id <- 1
cfg <- config::get(file = "omop/config.yml")
cd <- DatabaseConnector::createConnectionDetails(
  dbms = cfg$dbms, server = cfg$server, port = cfg$port,
  user = cfg$user, password = cfg$password,
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)
con <- DatabaseConnector::connect(cd)

# fill in the schemas, then turn the OHDSI SQL into PostgreSQL
sql <- SqlRender::readSql("omop/cohort/nsclc_platinum.sql")
sql <- SqlRender::render(sql, cdm_schema = cfg$cdm_schema,
                         results_schema = cfg$results_schema, cohort_id = cohort_id)
sql <- SqlRender::translate(sql, targetDialect = cfg$dbms)
DatabaseConnector::executeSql(con, sql, progressBar = FALSE)

q <- function(sql) {
  sql <- SqlRender::render(sql, cdm = cfg$cdm_schema, res = cfg$results_schema, id = cohort_id,
                         warnOnMissingParameters = FALSE)
  as_tibble(DatabaseConnector::querySql(con, SqlRender::translate(sql, cfg$dbms), snakeCaseToCamelCase = TRUE))
}

# ---- attrition: rules applied one after the other, so each patient is
# counted at the first rule they fail ----
cand <- q("SELECT * FROM @res.nsclc_candidates")
rules <- c(nsclc180d      = "NSCLC diagnosed in the 180 days before the first platinum dose",
           noSclc         = "Never diagnosed with small cell lung cancer",
           paclitaxel7d   = "Paclitaxel within 7 days of the first platinum dose",
           noOtherCancer5y = "No other cancer in the previous 5 years",
           adult          = "Aged 18 or over at the first dose",
           priorObs365d   = "At least a year of records before the first dose")

keep <- rep(TRUE, nrow(cand))
attrition <- tibble(step = "First platinum dose (cisplatin or carboplatin)", patients = nrow(cand))
for (r in names(rules)) {
  keep <- keep & cand[[r]] == 1
  attrition <- add_row(attrition, step = rules[[r]], patients = sum(keep))
}
attrition <- mutate(attrition, removed = lag(patients) - patients)
print(attrition)
write.csv(attrition, "omop/cohort/attrition.csv", row.names = FALSE)

# ---- who is in the cohort ----
# follow-up runs from the first dose to death, or to the end of the
# patient's records (censored); deaths after that date can't be seen anyway
coh <- q("
  SELECT c.subject_id, c.cohort_start_date, c.cohort_end_date, cand.age,
         p.gender_source_value AS sex, d.death_date
  FROM @res.cohort c
  JOIN @res.nsclc_candidates cand ON cand.person_id = c.subject_id
  JOIN @cdm.person p ON p.person_id = c.subject_id
  LEFT JOIN @cdm.death d ON d.person_id = c.subject_id
  WHERE c.cohort_definition_id = @id") %>%
  mutate(died = !is.na(deathDate) & deathDate <= cohortEndDate,
         end = if_else(died, deathDate, cohortEndDate),
         months = as.numeric(end - cohortStartDate) / 30.4375)
DatabaseConnector::disconnect(con)

stopifnot(nrow(coh) == tail(attrition$patients, 1), !anyDuplicated(coh$subjectId))

summ <- tibble(
  measure = c("Patients", "Age at first dose, median (IQR)", "Age 65 or over, %", "Male, %",
              "First dose, years", "Deaths", "Follow-up in months, median (range)"),
  value = c(nrow(coh),
            sprintf("%.0f (%.0f-%.0f)", median(coh$age), quantile(coh$age, .25), quantile(coh$age, .75)),
            sprintf("%.1f", 100 * mean(coh$age >= 65)),
            sprintf("%.1f", 100 * mean(coh$sex == "M")),
            paste(range(format(coh$cohortStartDate, "%Y")), collapse = "-"),
            sprintf("%d (%.1f%%)", sum(coh$died), 100 * mean(coh$died)),
            sprintf("%.1f (%.1f-%.1f)", median(coh$months), min(coh$months), max(coh$months))))
print(summ)
write.csv(summ, "omop/cohort/cohort_summary.csv", row.names = FALSE)
