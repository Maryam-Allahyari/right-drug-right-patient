# Who is in the cohort, described the standard OHDSI way: FeatureExtraction
# builds covariates for every patient from the OMOP tables (demographics,
# conditions and drugs in the year before the first platinum dose, the
# Charlson comorbidity index) and summarises them as a "Table 1".
# The same code works on any OMOP database, so a real-world cohort would be
# described exactly like this one.
#
# Before running:
#   - the cohort built with omop/cohort/build_cohort.R
#   - FeatureExtraction, DatabaseConnector, dplyr and config installed
#
# From the repo root:
#   Rscript omop/analysis/characterise_cohort.R
#
# Output in omop/results/: table1_cohort.csv and top_conditions_prior_year.csv.

suppressPackageStartupMessages(library(dplyr))

out <- "omop/results"
cfg <- config::get(file = "omop/config.yml")
cd <- DatabaseConnector::createConnectionDetails(
  dbms = cfg$dbms, server = cfg$server, port = cfg$port,
  user = cfg$user, password = cfg$password,
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)

# everything is looked up in the year before the index date, index day included
settings <- FeatureExtraction::createCovariateSettings(
  useDemographicsGender = TRUE, useDemographicsAgeGroup = TRUE,
  useDemographicsIndexYear = TRUE,
  useConditionGroupEraLongTerm = TRUE, useDrugGroupEraLongTerm = TRUE,
  useCharlsonIndex = TRUE,
  longTermStartDays = -365, endDays = 0
)

# aggregated = one row per covariate (share of patients), not per patient
cov <- FeatureExtraction::getDbCovariateData(
  connectionDetails = cd, cdmDatabaseSchema = cfg$cdm_schema,
  cohortDatabaseSchema = cfg$results_schema, cohortTable = "cohort", cohortIds = 1,
  covariateSettings = settings, aggregated = TRUE
)

# ---- Table 1 with FeatureExtraction's standard layout ----
t1 <- FeatureExtraction::createTable1(cov, output = "one column")
print(t1)
write.csv(t1, file.path(out, "table1_cohort.csv"), row.names = FALSE)

# ---- most common conditions in the year before the first dose ----
# condition group eras roll each diagnosis up the SNOMED hierarchy, so broad
# groups ("Disorder of lung") sit next to specific ones. Concepts every
# patient has (the lung cancer itself and its parent groups) are left out
ref <- collect(cov$covariateRef)
top <- collect(cov$covariates) %>%
  inner_join(ref, by = "covariateId") %>%
  filter(analysisId == 210) %>%   # 210 = ConditionGroupEraLongTerm
  transmute(condition = sub(".*relative to index: ", "", covariateName),
            patients = sumValue, percent = round(100 * averageValue, 1)) %>%
  filter(percent < 100) %>%
  arrange(desc(patients)) %>%
  head(25)
print(top, n = 25)
write.csv(top, file.path(out, "top_conditions_prior_year.csv"), row.names = FALSE)
