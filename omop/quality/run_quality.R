# Data quality checks on the OMOP database before any cohort is built.
#   - Achilles: descriptive statistics for every table (how many people,
#     records per year, most common concepts...). The results go into the
#     "results" schema, where ATLAS and the DQD read them.
#   - DataQualityDashboard (DQD): about 2,300 standard checks (on this database) of the
#     conformance, completeness and plausibility of the data, each with a
#     threshold for the percentage of rows allowed to fail.
#
# Before running:
#   - the OMOP database loaded with omop/etl/load_omop.R
#   - Achilles, DataQualityDashboard and jsonlite installed
#
# From the repo root:
#   Rscript omop/quality/run_quality.R
#
# Takes about 15 minutes. The full output (Achilles logs and the DQD json)
# goes to M:/Data/omop/quality, the summaries in omop/quality/ are small
# enough to keep in git. To browse the checks:
#   DataQualityDashboard::viewDqDashboard("M:/Data/omop/quality/dqd_results.json")

library(dplyr)

out_dir <- "M:/Data/omop/quality"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- config::get(file = "omop/config.yml")
cd <- DatabaseConnector::createConnectionDetails(
  dbms = cfg$dbms, server = cfg$server, port = cfg$port,
  user = cfg$user, password = cfg$password,
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)
source_name <- "Synthea lung cancer population"

step <- function(what) cat("
==", what, format(Sys.time(), "%H:%M"), "
")

# Both tools run on one thread. With several, each ran into a different
# problem on this setup: Windows had reserved the port the workers wanted,
# the workers couldn't read the connection details, Achilles' workers wrote
# to its log at the same moment and broke a line it reads back, and the DQD
# tries to close a connection it never opened. On a database this size one
# thread takes only a few minutes more.
step("Achilles")
Achilles::achilles(cd, cdmDatabaseSchema = cfg$cdm_schema, resultsDatabaseSchema = cfg$results_schema,
                   sourceName = source_name, cdmVersion = cfg$cdm_version,
                   numThreads = 1, outputFolder = file.path(out_dir, "achilles"))

# The vocabulary tables are skipped (the DQD default), they come from Athena
# as they are. So are the cohort tables, which don't exist until a cohort is
# built; checking them only gives SQL errors.
step("Data Quality Dashboard")
skip_tables <- c("CONCEPT", "VOCABULARY", "CONCEPT_ANCESTOR", "CONCEPT_RELATIONSHIP",
                 "CONCEPT_CLASS", "CONCEPT_SYNONYM", "RELATIONSHIP", "DOMAIN", "PACK_CONTENT",
                 "CONCEPT_METADATA", "CONCEPT_RELATIONSHIP_METADATA",
                 "COHORT", "COHORT_DEFINITION")
DataQualityDashboard::executeDqChecks(
  cd, cdmDatabaseSchema = cfg$cdm_schema, resultsDatabaseSchema = cfg$results_schema,
  cdmSourceName = source_name, cdmVersion = cfg$cdm_version, numThreads = 1,
  outputFolder = out_dir, outputFile = "dqd_results.json",
  tablesToExclude = skip_tables
)

# ---- summaries for the repo ----
step("summaries")
res <- jsonlite::fromJSON(file.path(out_dir, "dqd_results.json"))
checks <- as_tibble(res$CheckResults)

# pass/fail per category, the same numbers as the dashboard's overview page.
# "Not applicable" means the table or field is empty, so there's nothing to check;
# errors are checks whose SQL didn't run
overview <- checks %>%
  mutate(status = case_when(isError == 1 ~ "error", notApplicable == 1 ~ "not applicable",
                            failed == 1 ~ "failed", TRUE ~ "passed")) %>%
  count(category, context, status) %>%
  tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)
print(overview)
write.csv(overview, "omop/quality/dqd_overview.csv", row.names = FALSE)

# every failed check, worst first
failed <- checks %>%
  filter(failed == 1) %>%
  select(category, subcategory, context, checkLevel, checkName, cdmTableName, cdmFieldName,
         conceptId, numViolatedRows, numDenominatorRows, pctViolatedRows, thresholdValue,
         checkDescription) %>%
  arrange(desc(pctViolatedRows))
write.csv(failed, "omop/quality/dqd_failed.csv", row.names = FALSE)
cat(nrow(failed), "of", nrow(checks), "checks failed\n")
step("done")
