# Synthea → OMOP CDM v5.4: mapping notes

The conversion itself is done with [ETL-Synthea](https://github.com/OHDSI/ETL-Synthea), which already has a reviewed mapping for every table. Instead of redesigning it, I checked it against my own scan of the data (`../whiterabbit/ScanReport.xlsx`) and wrote down what matters for this project: what each table the cohort uses is built from, things that aren't obvious from the ETL documentation, and where information is lost.

The ETL-Synthea pages below are the source of truth for the full column-by-column logic.

## Tables the cohort uses

| OMOP table | Built from (Synthea) | Key logic | ETL-Synthea page |
|---|---|---|---|
| `person` | `patients.csv` | One row per patient; year, month and day of birth split from `BIRTHDATE`; gender, race and ethnicity mapped to standard concepts | [person](https://ohdsi.github.io/ETL-Synthea/articles/person.html) |
| `death` | `encounters.csv` | Encounters with code 308646001 (*Death certification*); death date = start of that encounter; cause of death from its reason code | [death](https://ohdsi.github.io/ETL-Synthea/articles/death.html) |
| `observation_period` | `encounters.csv` | One period per patient, from their first encounter start to their last encounter end | [observation_period](https://ohdsi.github.io/ETL-Synthea/articles/observation_period.html) |
| `visit_occurrence` | `encounters.csv` | Encounter class → inpatient (9201), outpatient (9202), emergency (9203); inpatient stays less than a day apart are merged into one visit | [visit_occurrence](https://ohdsi.github.io/ETL-Synthea/articles/visit_occurrence.html) |
| `condition_occurrence` | `conditions.csv` | SNOMED code → standard SNOMED condition concept; start and stop dates copied | [condition_occurrence](https://ohdsi.github.io/ETL-Synthea/articles/condition_occurrence.html) |
| `drug_exposure` | `medications.csv`, `immunizations.csv` | RxNorm (medications) and CVX (vaccines) → standard drug concepts; if there's no stop date, end date = start date | [drug_exposure](https://ohdsi.github.io/ETL-Synthea/articles/drug_exposure.html) |
| `procedure_occurrence` | `procedures.csv` | SNOMED code → standard SNOMED procedure concept | [procedure_occurrence](https://ohdsi.github.io/ETL-Synthea/articles/procedure_occurrence.html) |

## What I found in the data

**Death date is not `patients.DEATHDATE`.** ETL-Synthea takes death from the *Death certification* encounter, not from the `DEATHDATE` column. In this data every patient with a `DEATHDATE` (748) has exactly one such encounter, but it is dated 0 to 14 days *after* `DEATHDATE` (median 4 days). So death dates in OMOP are a few days later than in the source. For survival measured in months this hardly matters, but it's the kind of thing to know before comparing with another source.

**Some conditions, including other cancers, are coded in ICD-10.** About 3% of condition records (1,644; 238 different codes) use ICD-10 instead of SNOMED:

| ICD-10 chapter | Records | Patients | What it is |
|---|---|---|---|
| K | 1,078 | 413 | Dental (caries, abscesses, tooth disorders) |
| C | 252 | 174 | Malignant tumours, mostly mouth, tongue, throat and larynx; 2 records of lung cancer (C34.9) |
| Z, M, D, R, S | 314 | | Check-ups, jaw disorders, benign tumours, injuries |

All 835 lung cancer diagnoses used for the cohort are SNOMED. But 142 of the NSCLC patients also have an ICD-10 cancer code, and the cohort excludes patients with another cancer in the previous five years. The ETL-Synthea documentation describes condition codes being looked up as SNOMED; if ICD-10 codes aren't mapped, these cancers end up with `condition_concept_id = 0` and the exclusion would silently miss them. **To check after loading:** how many condition records have concept 0, and whether the ICD-10 cancer codes are mapped to standard SNOMED concepts. If not, they have to be mapped through the vocabulary's ICD-10 → SNOMED *Maps to* relationships before building the cohort.

**Two date formats.** `patients.csv` and `conditions.csv` use plain dates (`2021-10-05`); `encounters`, `medications` and `procedures` use timestamps with a time zone (`2021-10-05T03:30:19Z`). OMOP keeps the date and sets times to midnight, so the time of day is lost. Not a problem for this project.

**Chemotherapy is recorded over and over.** `medications.csv` has 428,288 rows, mostly cisplatin and paclitaxel given repeatedly for years. Each row becomes a `drug_exposure` record, so that table will be big. 3,700 medication rows have no stop date; ETL-Synthea then sets the end date to the start date.

**Direct identifiers are dropped.** Names, SSN, driver's licence, passport, address, coordinates and income in `patients.csv` are not carried into OMOP. `person_source_value` keeps the Synthea patient ID, which is a random UUID.

## Where information is lost

- **Race:** only *white* (8527), *black* (8516) and *asian* (8515) get a concept. *Native* and *other* become 0, so about 7% of patients (62 of 835) have no race concept. The original value stays in `race_source_value`.
- **Ethnicity:** *hispanic* → 38003563; everything else (including *nonhispanic*) → 0. So "not Hispanic" and "unknown" look the same in `ethnicity_concept_id`. Again, `ethnicity_source_value` keeps the original.
- **Dental procedures:** CDT codes (about 3% of procedures) aren't in the Athena vocabularies, so these get `procedure_concept_id = 0`. The code itself stays in `procedure_source_value`.
- **Time of day** for encounters, drugs and procedures (see above).
- **Location, provider and care site** are mostly left empty for `person`; this project doesn't need them.

## What this means for the cohort

- **Start of follow-up** (first chemotherapy after diagnosis) comes from `drug_exposure` with the RxNorm codes for cisplatin, etoposide and paclitaxel.
- **Death** comes from the OMOP `death` table, i.e. the certification date (a few days after the true date).
- **End of follow-up** for patients still alive is the end of their `observation_period` (last encounter), which is what makes censoring possible.
- **SCLC** patients are excluded by their SNOMED code (254632001); NSCLC is 254637007.
- **Other cancers** used for the exclusion rule must include the ICD-10-coded ones (see above), not only SNOMED.
