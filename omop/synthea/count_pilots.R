# Quick look at the two Synthea pilots (lung and breast) to decide which
# cancer to build the cohort on. Reads the raw Synthea CSVs directly,
# no OMOP database needed at this stage.
#
# Before running:
#   - generate the pilots with run_pilots.sh (output goes to M:/Data/synthea)
#   - install.packages(c("DBI", "duckdb")) if you don't have them
#
# From the repo root:
#   Rscript omop/synthea/count_pilots.R
#
# Output: a small table with the lung and breast counts side by side.
# A pilot that hasn't been generated yet is just left out.

library(DBI)
library(duckdb)

base <- "M:/Data/synthea"   # change this if your Synthea output is somewhere else

# Codes that count as treatment. These are the ones the Synthea modules use,
# so if you switch Synthea version, check they haven't changed.
treat <- list(
  lung   = list(drugs = c("1736854", "226719", "583214"),     # cisplatin, etoposide, paclitaxel
                procs = c("703423002")),                       # combined chemo + radiation
  breast = list(drugs = c("1732186", "1790099", "583214"),     # epirubicin, doxorubicin, paclitaxel
                procs = c("367336001"))                        # chemotherapy
)

# Diagnosis codes (SNOMED). Lung has two types, we count them separately.
dx <- list(
  lung   = c(NSCLC = "254637007", SCLC = "254632001"),
  breast = c(breast = "254837009")
)

con <- dbConnect(duckdb(shared_home = FALSE))
# Some Synthea date columns are plain dates and some are timestamps with a
# time zone, and DuckDB won't cast the second kind straight to DATE.
# d() just takes the first 10 characters (yyyy-mm-dd).
invisible(dbExecute(con, "CREATE MACRO d(x) AS CAST(LEFT(CAST(x AS VARCHAR), 10) AS DATE)"))
q <- function(sql) dbGetQuery(con, sql)
in_list <- function(x) paste0("'", x, "'", collapse = ", ")   # c("a","b") -> 'a', 'b' for SQL

count_pilot <- function(name) {
  dir <- file.path(base, paste0("pilot_", name), "csv")
  # point DuckDB at the CSVs (nothing is copied, it reads them in place)
  for (t in c("patients", "conditions", "medications", "procedures")) {
    dbExecute(con, sprintf("CREATE OR REPLACE VIEW %s AS SELECT * FROM read_csv_auto('%s/%s.csv')",
                           t, dir, t))
  }

  # earliest diagnosis date per patient
  dbExecute(con, sprintf("
    CREATE OR REPLACE TEMP TABLE dx AS
    SELECT PATIENT, d(MIN(START)) AS dx_date, ANY_VALUE(CODE) AS dx_code
    FROM conditions WHERE CAST(CODE AS VARCHAR) IN (%s) GROUP BY PATIENT", in_list(dx[[name]])))

  # treatment records from the diagnosis date onwards
  # (anything before diagnosis can't be treatment for this cancer)
  dbExecute(con, sprintf("
    CREATE OR REPLACE TEMP TABLE tx AS
    SELECT m.PATIENT, d(m.START) AS tx_date, CAST(m.CODE AS VARCHAR) AS code, 'drug' AS kind
    FROM medications m JOIN dx USING (PATIENT)
    WHERE CAST(m.CODE AS VARCHAR) IN (%s) AND d(m.START) >= dx.dx_date
    UNION ALL
    SELECT p.PATIENT, d(p.START), CAST(p.CODE AS VARCHAR), 'procedure'
    FROM procedures p JOIN dx USING (PATIENT)
    WHERE CAST(p.CODE AS VARCHAR) IN (%s) AND d(p.START) >= dx.dx_date",
    in_list(treat[[name]]$drugs), in_list(treat[[name]]$procs)))

  # second course: a gap of more than 90 days between treatment records.
  # 90 days is a rough rule, Synthea doesn't record treatment lines.
  courses <- q("
    SELECT COUNT(*) AS n FROM (
      SELECT PATIENT, SUM(new_course) AS courses FROM (
        SELECT PATIENT, CASE WHEN tx_date - LAG(tx_date) OVER (PARTITION BY PATIENT ORDER BY tx_date) > 90
                             THEN 1 ELSE 0 END AS new_course
        FROM (SELECT DISTINCT PATIENT, tx_date FROM tx)) GROUP BY PATIENT)
    WHERE courses >= 1")$n

  # deaths within 30 months of the first treatment.
  # 913 days = 30 months, same follow-up cut-off as the simulated trial.
  first_tx <- "SELECT PATIENT, MIN(tx_date) AS tx1 FROM tx GROUP BY PATIENT"
  deaths30 <- q(sprintf("
    SELECT COUNT(*) AS n FROM (%s) f JOIN patients p ON p.Id = f.PATIENT
    WHERE p.DEATHDATE IS NOT NULL
      AND d(p.DEATHDATE) - f.tx1 BETWEEN 0 AND 913", first_tx))$n

  # one row of counts for this pilot
  data.frame(
    pilot = name,
    patients = q("SELECT COUNT(*) AS n FROM patients")$n,
    dead = q("SELECT COUNT(*) AS n FROM patients WHERE DEATHDATE IS NOT NULL")$n,
    with_diagnosis = q("SELECT COUNT(*) AS n FROM dx")$n,
    diagnosis_split = paste(apply(q("SELECT dx_code, COUNT(*) AS n FROM dx GROUP BY dx_code"), 1,
                                  paste, collapse = "="), collapse = "; "),
    treated_after_dx = q("SELECT COUNT(DISTINCT PATIENT) AS n FROM tx")$n,
    treatment_records = q("SELECT COUNT(*) AS n FROM tx")$n,
    distinct_drugs = q("SELECT COUNT(DISTINCT code) AS n FROM tx WHERE kind = 'drug'")$n,
    deaths_30m_after_tx1 = deaths30,
    second_course = courses
  )
}

res <- do.call(rbind, lapply(c("lung", "breast"), function(n) {
  if (dir.exists(file.path(base, paste0("pilot_", n), "csv"))) count_pilot(n) else NULL
}))
print(t(res))   # transposed so the two pilots sit side by side
dbDisconnect(con, shutdown = TRUE)
