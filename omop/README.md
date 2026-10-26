# Real-world cohort in the OMOP CDM

This folder holds the real-world part of the project: synthetic patient records (Synthea) converted to the OMOP Common Data Model, a cohort of treated cancer patients, and its survival compared with the simulated trial's control arm.

## Choosing the disease (pilot, 2026-10-26)

The analysis needs three things from the data: chemotherapy records, to define when follow-up starts; deaths, for overall survival; and enough patients in a general (not restricted) population. Synthea 3.3.0 has six cancer-related modules:

| Module | Decision | Reason |
|---|---|---|
| Lung cancer | Piloted | Chemotherapy, staging and deaths are modelled |
| Breast cancer | Piloted | Chemotherapy and deaths are modelled |
| Colorectal cancer | Not piloted in this round | Has chemotherapy (oxaliplatin, leucovorin) and deaths; could be compared later |
| Acute myeloid leukaemia | Excluded | No chemotherapy is modelled, so there is no treatment start to use as time zero |
| Veteran lung cancer | Excluded | Same model as lung cancer, but only for the simulated veteran population |
| Veteran prostate cancer | Excluded | Only for the simulated veteran population, and men only |

Two small Synthea 3.3.0 populations were generated, keeping only patients with the cancer diagnosis (`synthea/keep_lung.json`, `synthea/keep_breast.json`; commands in `synthea/run_pilots.sh`; seed 2026, ages 45–90). Counts were taken straight from the Synthea CSVs with `synthea/count_pilots.R`.

| | Lung | Breast |
|---|---|---|
| Patients | 204 (NSCLC 169, SCLC 35) | 134 |
| Treated with chemotherapy after diagnosis | 204 | 41 |
| Deaths within 30 months of first treatment | 39 | 9 |
| Patients with a second treatment course (> 90-day gap) | 14 | 0 |
| Drugs | cisplatin, etoposide, paclitaxel | epirubicin, doxorubicin, paclitaxel |

**Decision: lung cancer, non-small cell (NSCLC) only.** Every lung patient has treatment records, there are about four times as many deaths within 30 months, and only lung has second treatment courses. Small cell lung cancer is a different disease and is left out.

**Things to keep in mind**
- In Synthea, lung cancer patients live much longer than in reality: the median time from first treatment to death is about 45 months, while real advanced NSCLC has a median survival of about a year. Chemotherapy records also continue for years. These are features of the simulator, not findings.
- Synthea generates men and women equally, so many attempts to find a breast cancer patient failed (the "failed to produce a matching patient" messages in the breast log).
- `-p` in Synthea counts living patients only; deceased patients come on top.

## Full population

`synthea/run_full.sh` generates the lung cancer population with the same settings as the pilot (seed 2026, ages 45–90), asking for 120 living patients:

| | NSCLC | SCLC |
|---|---|---|
| Patients | 714 | 121 |
| Treated with chemotherapy after diagnosis | 713 | 121 |
| Deaths within 30 months of first treatment | 123 | 22 |

SCLC patients stay in the data and are excluded when the cohort is built.

## Source data profile (WhiteRabbit)

Before converting anything to OMOP, the Synthea CSVs were scanned with WhiteRabbit v1.0.0 (settings in `whiterabbit/scan.ini`). The scan report lists every table and column with its type, how often it is empty and its most common values. Values seen fewer than 5 times are left out of the report, and for large tables the frequencies come from the first 100,000 rows only.

What matters for the conversion:
- **Code systems are mixed.** Conditions are mostly SNOMED CT, but about 3% are ICD-10; procedures are SNOMED CT plus about 3% CDT (dental) codes; medications are RxNorm. ICD-10 codes have to be mapped to standard SNOMED concepts.
- **Two date formats.** Patients and conditions use plain dates (`2021-10-05`); encounters, medications and procedures use timestamps with a time zone (`2021-10-05T03:30:19Z`).
- **Death is a column, not a table.** `DEATHDATE` is in `patients.csv` (filled for 748 of 835 patients) and becomes the OMOP `death` table.
- **Personal identifiers** (names, SSN, address, coordinates; all synthetic) are in `patients.csv`. They are not carried into OMOP.
- **Chemotherapy dominates the medication records:** cisplatin and paclitaxel are about 73% of the first 100,000 rows, because Synthea records them repeatedly for years. The OMOP `drug_exposure` table will be large.
- About 29% of conditions have no end date (ongoing conditions).

## Data location
Synthea output is written to `M:/Data/synthea/` and is not tracked in git; it can be regenerated with the commands in `synthea/`.
