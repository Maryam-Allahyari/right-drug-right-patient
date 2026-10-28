# Loads the full lung cancer Synthea population into an OMOP CDM v5.4
# database with ETL-Synthea: raw Synthea tables go into the "native" schema,
# the OMOP tables and the vocabulary into "cdm".
#
# Before running:
#   - PostgreSQL with the synthea_omop database (omop/sql/00_setup_database.sql)
#   - omop/config.yml with your password (copy config.yml.example)
#   - the Athena vocabulary unzipped in M:/Data/vocab/v5 (list in omop/README.md)
#   - the Synthea CSVs from synthea/run_full.sh
#   - ETLSyntheaBuilder, DatabaseConnector and config installed, and the
#     PostgreSQL JDBC driver in the folder DATABASECONNECTOR_JAR_FOLDER points to
#
# From the repo root:
#   Rscript omop/etl/load_omop.R
#
# Starts from empty cdm and native schemas every time, so it can just be
# rerun. Takes about 35 minutes, most of it loading the vocabulary.
# Close Docker or other memory-hungry programs first, the vocabulary load
# reads whole files into memory.

library(ETLSyntheaBuilder)

# where things are on my machine, change these to yours
synthea_csv <- "M:/Data/synthea/full_lung/csv"
vocab_csv   <- "M:/Data/vocab/v5"
# bulk loading hands the data to psql instead of sending it row by row
# through Java, which is the only way the vocabulary loads in reasonable time
Sys.setenv(POSTGRES_PATH = "C:/Program Files/PostgreSQL/18/bin")

cfg <- config::get(file = "omop/config.yml")
cd <- DatabaseConnector::createConnectionDetails(
  dbms = cfg$dbms, server = cfg$server, port = cfg$port,
  user = cfg$user, password = cfg$password,
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER")
)

step <- function(what) cat("\n==", what, format(Sys.time(), "%H:%M"), "\n")

# empty the two schemas first, the create functions fail if the tables exist
step("reset schemas")
con <- DatabaseConnector::connect(cd)
for (s in c(cfg$cdm_schema, cfg$synthea_schema)) {
  DatabaseConnector::executeSql(con, sprintf("DROP SCHEMA IF EXISTS %s CASCADE; CREATE SCHEMA %s;", s, s))
}
DatabaseConnector::disconnect(con)

step("create tables")
CreateCDMTables(cd, cdmSchema = cfg$cdm_schema, cdmVersion = cfg$cdm_version)
CreateSyntheaTables(cd, syntheaSchema = cfg$synthea_schema, syntheaVersion = cfg$synthea_version)

# LoadSyntheaTables loads every CSV in the folder, so the ones the ETL doesn't
# need are left out: imaging_studies isn't used at all and claims_transactions
# only feeds the cost table. Together that's 5 of the 6 GB.
# Their native tables stay empty, which the ETL is fine with.
step("load Synthea CSVs")
skip <- c("imaging_studies.csv", "claims_transactions.csv")
load_dir <- file.path(dirname(synthea_csv), "csv_for_etl")
dir.create(load_dir, showWarnings = FALSE)
keep <- setdiff(list.files(synthea_csv, pattern = "\\.csv$"), skip)
file.copy(file.path(synthea_csv, keep), load_dir, overwrite = TRUE)
LoadSyntheaTables(cd, syntheaSchema = cfg$synthea_schema, syntheaFileLoc = load_dir, bulkLoad = TRUE)

step("load vocabulary")
LoadVocabFromCsv(cd, cdmSchema = cfg$cdm_schema, vocabFileLoc = vocab_csv, bulkLoad = TRUE)

# lookup tables from source codes to standard concepts, plus visit rollup;
# the extra indices make the next step a lot faster
step("map and rollup tables")
CreateMapAndRollupTables(cd, cdmSchema = cfg$cdm_schema, syntheaSchema = cfg$synthea_schema,
                         cdmVersion = cfg$cdm_version, syntheaVersion = cfg$synthea_version)
CreateExtraIndices(cd, cdmSchema = cfg$cdm_schema, syntheaSchema = cfg$synthea_schema,
                   syntheaVersion = cfg$synthea_version)

# With PostgreSQL's defaults the observation step ran for 40 minutes without
# writing a row: the planner expected about 20 matching rows (there are
# 470,000), picked row-by-row lookups and spilled to disk with 4 MB of memory.
# More memory and hash joins only, just while the ETL runs: 8 seconds.
# Set on the database because LoadEventTables opens its own connection.
db <- sub(".*/", "", cfg$server)
set_db <- function(sql) {
  con <- DatabaseConnector::connect(cd)
  DatabaseConnector::executeSql(con, sprintf(sql, db, db), progressBar = FALSE)
  DatabaseConnector::disconnect(con)
}
set_db("ALTER DATABASE %s SET work_mem = '1GB'; ALTER DATABASE %s SET enable_nestloop = off;")

step("fill OMOP tables")
LoadEventTables(cd, cdmSchema = cfg$cdm_schema, syntheaSchema = cfg$synthea_schema,
                cdmVersion = cfg$cdm_version, syntheaVersion = cfg$synthea_version,
                cdmSourceName = "Synthea lung cancer population",
                cdmSourceAbbreviation = "Synthea lung",
                cdmHolder = "right-drug-right-patient",
                createIndices = TRUE)
set_db("ALTER DATABASE %s RESET work_mem; ALTER DATABASE %s RESET enable_nestloop;")   # back to normal

# quick sanity check: rows in the tables the cohort uses
step("row counts")
con <- DatabaseConnector::connect(cd)
tabs <- c("person", "observation_period", "visit_occurrence", "condition_occurrence",
          "drug_exposure", "procedure_occurrence", "death", "concept")
counts <- sapply(tabs, function(t)
  DatabaseConnector::querySql(con, sprintf("SELECT COUNT(*) AS n FROM %s.%s", cfg$cdm_schema, t))[[1]])
print(data.frame(rows = counts))
DatabaseConnector::disconnect(con)
step("done")
