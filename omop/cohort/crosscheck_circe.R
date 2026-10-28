# The same cohort built a second way, with the standard OHDSI cohort tools,
# as a check on the hand-written SQL (nsclc_platinum.sql):
#   - Capr writes the definition in R and exports it as Circe JSON, the format
#     ATLAS uses (nsclc_platinum.json can be imported into ATLAS as it is)
#   - CirceR turns the JSON into SQL
#   - CohortGenerator runs it and counts patients per inclusion rule
# Then the two cohorts are compared patient by patient.
#
# Before running:
#   - omop/cohort/build_cohort.R run first (the SQL cohort)
#   - Capr (from github.com/OHDSI/Capr), CirceR, CohortGenerator,
#     DatabaseConnector, dplyr and config installed
#
# From the repo root:
#   Rscript omop/cohort/crosscheck_circe.R
#
# Output: nsclc_platinum.json, attrition_circe.csv and crosscheck.csv in
# omop/cohort/, and the cohort in results.cohort_circe.

suppressPackageStartupMessages({library(dplyr); library(Capr)})

cfg <- config::get(file = "omop/config.yml")
cd <- DatabaseConnector::createConnectionDetails(
  dbms = cfg$dbms, server = cfg$server, port = cfg$port,
  user = cfg$user, password = cfg$password,
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)

# ---- concept sets: the same concept IDs as in the SQL, with descendants ----
platinum   <- cs(descendants(1397599, 1344905), name = "Platinum (cisplatin, carboplatin)")
paclitaxel <- cs(descendants(1378382), name = "Paclitaxel")
nsclc      <- cs(descendants(4115276), name = "Non-small cell lung cancer")
sclc       <- cs(descendants(4110591), name = "Small cell lung cancer")
other_cancer <- cs(descendants(443392), exclude(descendants(443388, 4112752, 4111921)),
                   name = "Other cancer (not lung, basal or squamous cell skin)")

# ---- the definition ----
# Circe checks prior observation at entry, so the "a year of records before"
# rule sits in the entry event here rather than in the list of rules.
# Age is checked on the entry drug record itself (age at the first dose).
rules <- list(
  "NSCLC diagnosed in the 180 days before the first platinum dose" =
    withAll(atLeast(1, conditionOccurrence(nsclc), duringInterval(eventStarts(-180, 0)))),
  "Never diagnosed with small cell lung cancer" =
    withAll(exactly(0, conditionOccurrence(sclc),
                    duringInterval(eventStarts(-Inf, Inf), ignoreObservationPeriod = TRUE))),
  "Paclitaxel within 7 days of the first platinum dose" =
    withAll(atLeast(1, drugExposure(paclitaxel), duringInterval(eventStarts(-7, 7)))),
  "No other cancer in the previous 5 years" =
    withAll(exactly(0, conditionOccurrence(other_cancer), duringInterval(eventStarts(-1825, 0)))),
  "Aged 18 or over at the first dose" =
    withAll(atLeast(1, drugExposure(platinum, age(gte(18))), duringInterval(eventStarts(0, 0))))
)
def <- cohort(
  entry = entry(drugExposure(platinum, firstOccurrence()),
                observationWindow = continuousObservation(priorDays = 365),
                primaryCriteriaLimit = "First"),
  attrition = do.call(attrition, rules),
  exit = exit(endStrategy = observationExit())
)
json <- as.json(def)
writeLines(json, "omop/cohort/nsclc_platinum.json")

# ---- generate it with CirceR and CohortGenerator ----
cohort_set <- data.frame(
  cohortId = 2, cohortName = "NSCLC first platinum (Circe)", json = json,
  sql = CirceR::buildCohortQuery(CirceR::cohortExpressionFromJson(json),
                                 CirceR::createGenerateOptions(generateStats = TRUE))
)
tables <- CohortGenerator::getCohortTableNames(cohortTable = "cohort_circe")
CohortGenerator::createCohortTables(cd, cohortDatabaseSchema = cfg$results_schema,
                                    cohortTableNames = tables)
CohortGenerator::generateCohortSet(cd, cdmDatabaseSchema = cfg$cdm_schema,
                                   cohortDatabaseSchema = cfg$results_schema,
                                   cohortTableNames = tables, cohortDefinitionSet = cohort_set,
                                   incremental = FALSE)

# Attrition. Circe stores, for every combination of passed rules, how many
# people have it: rule i is bit i of inclusionRuleMask (31 = all five passed).
# Counting the people who passed rules 1..k gives the same step-by-step
# attrition as build_cohort.R. Mode 0 = per person.
stats <- CohortGenerator::getCohortStats(cd, databaseId = "synthea",
                                         cohortDatabaseSchema = cfg$results_schema,
                                         cohortTableNames = tables)
masks <- filter(stats$cohortInclusionResultTable, modeId == 0)
passed_first_k <- sapply(seq_along(rules), function(k) {
  m <- 2^k - 1
  sum(masks$personCount[bitwAnd(masks$inclusionRuleMask, m) == m])
})
entry_n <- filter(stats$cohortSummaryStatsTable, modeId == 0)$baseCount
attrition_circe <- tibble(step = c("First platinum dose, with a year of records before", names(rules)),
                          patients = c(entry_n, passed_first_k)) %>%
  mutate(removed = lag(patients) - patients)
print(attrition_circe)
write.csv(attrition_circe, "omop/cohort/attrition_circe.csv", row.names = FALSE)

# ---- compare the two cohorts patient by patient ----
con <- DatabaseConnector::connect(cd)
both <- DatabaseConnector::querySql(con, SqlRender::translate(SqlRender::render("
  SELECT COALESCE(s.subject_id, c.subject_id) AS subject_id,
         s.cohort_start_date AS sql_start, c.cohort_start_date AS circe_start,
         s.cohort_end_date AS sql_end, c.cohort_end_date AS circe_end
  FROM (SELECT * FROM @res.cohort WHERE cohort_definition_id = 1) s
  FULL OUTER JOIN (SELECT * FROM @res.cohort_circe WHERE cohort_definition_id = 2) c
    ON c.subject_id = s.subject_id", res = cfg$results_schema), cfg$dbms),
  snakeCaseToCamelCase = TRUE)
DatabaseConnector::disconnect(con)

check <- tibble(
  measure = c("Patients in the SQL cohort", "Patients in the Circe cohort", "In both",
              "Only in the SQL cohort", "Only in the Circe cohort",
              "In both, different start date", "In both, different end date"),
  n = c(sum(!is.na(both$sqlStart)), sum(!is.na(both$circeStart)),
        sum(!is.na(both$sqlStart) & !is.na(both$circeStart)),
        sum(!is.na(both$sqlStart) & is.na(both$circeStart)),
        sum(is.na(both$sqlStart) & !is.na(both$circeStart)),
        sum(both$sqlStart != both$circeStart, na.rm = TRUE),
        sum(both$sqlEnd != both$circeEnd, na.rm = TRUE)))
print(check)
write.csv(check, "omop/cohort/crosscheck.csv", row.names = FALSE)
