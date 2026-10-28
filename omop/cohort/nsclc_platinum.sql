/*
NSCLC patients starting platinum chemotherapy for the first time:
the routine-care counterpart of the trial's control arm.

Written in OHDSI SQL (SQL Server style) and translated with SqlRender, so it
runs on any database the OMOP tools support. Parameters:
  @cdm_schema      OMOP tables and vocabulary
  @results_schema  where the candidates and the cohort table go
  @cohort_id       cohort_definition_id in the cohort table

Index date: first ever platinum dose (cisplatin or carboplatin).
Every person with a platinum dose is a candidate; each inclusion rule is a
0/1 column, so the R script can count how many drop out at each rule.
Concept IDs are standard concepts, always used with all their descendants.
*/

-- ---- concept sets ----
DROP TABLE IF EXISTS #codesets;
CREATE TABLE #codesets (codeset VARCHAR(20), concept_id BIGINT);

INSERT INTO #codesets
-- 1397599 cisplatin, 1344905 carboplatin (ingredients)
SELECT 'platinum', descendant_concept_id FROM @cdm_schema.concept_ancestor WHERE ancestor_concept_id IN (1397599, 1344905)
UNION ALL
-- 1378382 paclitaxel, the trial's other control drug
SELECT 'paclitaxel', descendant_concept_id FROM @cdm_schema.concept_ancestor WHERE ancestor_concept_id = 1378382
UNION ALL
-- 4115276 Non-small cell lung cancer (SNOMED 254637007)
SELECT 'nsclc', descendant_concept_id FROM @cdm_schema.concept_ancestor WHERE ancestor_concept_id = 4115276
UNION ALL
-- 4110591 Small cell carcinoma of lung (SNOMED 254632001)
SELECT 'sclc', descendant_concept_id FROM @cdm_schema.concept_ancestor WHERE ancestor_concept_id = 4110591
UNION ALL
-- any other cancer: 443392 Malignant neoplastic disease, minus lung cancer itself
-- (443388) and the skin cancers trials usually allow (4112752 basal cell,
-- 4111921 squamous cell carcinoma of skin)
SELECT 'other_cancer', descendant_concept_id FROM @cdm_schema.concept_ancestor
WHERE ancestor_concept_id = 443392
  AND descendant_concept_id NOT IN (
    SELECT descendant_concept_id FROM @cdm_schema.concept_ancestor
    WHERE ancestor_concept_id IN (443388, 4112752, 4111921));

-- ---- entry event: first platinum dose, inside an observation period ----
DROP TABLE IF EXISTS #first_platinum;
SELECT de.person_id, MIN(de.drug_exposure_start_date) AS index_date
INTO #first_platinum
FROM @cdm_schema.drug_exposure de
JOIN #codesets cs ON cs.concept_id = de.drug_concept_id AND cs.codeset = 'platinum'
GROUP BY de.person_id;

-- ---- one row per candidate with a flag per inclusion rule ----
DROP TABLE IF EXISTS @results_schema.nsclc_candidates;
SELECT fp.person_id, fp.index_date,
  op.observation_period_start_date AS obs_start,
  op.observation_period_end_date AS obs_end,
  YEAR(fp.index_date) - p.year_of_birth AS age,

  -- NSCLC diagnosed in the 180 days up to the first dose, so the platinum
  -- is the first treatment for this cancer (first line)
  CASE WHEN EXISTS (
    SELECT 1 FROM @cdm_schema.condition_occurrence co
    JOIN #codesets cs ON cs.concept_id = co.condition_concept_id AND cs.codeset = 'nsclc'
    WHERE co.person_id = fp.person_id
      AND co.condition_start_date BETWEEN DATEADD(day, -180, fp.index_date) AND fp.index_date
  ) THEN 1 ELSE 0 END AS nsclc_180d,

  -- never diagnosed with small cell lung cancer, a different disease
  CASE WHEN NOT EXISTS (
    SELECT 1 FROM @cdm_schema.condition_occurrence co
    JOIN #codesets cs ON cs.concept_id = co.condition_concept_id AND cs.codeset = 'sclc'
    WHERE co.person_id = fp.person_id
  ) THEN 1 ELSE 0 END AS no_sclc,

  -- paclitaxel within a week of the first platinum dose: the trial's control
  -- arm is cisplatin + paclitaxel
  CASE WHEN EXISTS (
    SELECT 1 FROM @cdm_schema.drug_exposure de
    JOIN #codesets cs ON cs.concept_id = de.drug_concept_id AND cs.codeset = 'paclitaxel'
    WHERE de.person_id = fp.person_id
      AND de.drug_exposure_start_date BETWEEN DATEADD(day, -7, fp.index_date) AND DATEADD(day, 7, fp.index_date)
  ) THEN 1 ELSE 0 END AS paclitaxel_7d,

  -- no other cancer in the 5 years before, the usual trial exclusion
  CASE WHEN NOT EXISTS (
    SELECT 1 FROM @cdm_schema.condition_occurrence co
    JOIN #codesets cs ON cs.concept_id = co.condition_concept_id AND cs.codeset = 'other_cancer'
    WHERE co.person_id = fp.person_id
      AND co.condition_start_date BETWEEN DATEADD(day, -1825, fp.index_date) AND fp.index_date
  ) THEN 1 ELSE 0 END AS no_other_cancer_5y,

  CASE WHEN YEAR(fp.index_date) - p.year_of_birth >= 18 THEN 1 ELSE 0 END AS adult,

  -- a year of records before the index date, so the look-back rules can see something
  CASE WHEN DATEDIFF(day, op.observation_period_start_date, fp.index_date) >= 365
       THEN 1 ELSE 0 END AS prior_obs_365d
INTO @results_schema.nsclc_candidates
FROM #first_platinum fp
JOIN @cdm_schema.person p ON p.person_id = fp.person_id
JOIN @cdm_schema.observation_period op
  ON op.person_id = fp.person_id
 AND fp.index_date BETWEEN op.observation_period_start_date AND op.observation_period_end_date;

-- ---- the cohort: candidates passing every rule, followed to the end of their records ----
IF OBJECT_ID('@results_schema.cohort', 'U') IS NULL
CREATE TABLE @results_schema.cohort (
  cohort_definition_id BIGINT, subject_id BIGINT,
  cohort_start_date DATE, cohort_end_date DATE);

DELETE FROM @results_schema.cohort WHERE cohort_definition_id = @cohort_id;

INSERT INTO @results_schema.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
SELECT @cohort_id, person_id, index_date, obs_end
FROM @results_schema.nsclc_candidates
WHERE nsclc_180d = 1 AND no_sclc = 1 AND paclitaxel_7d = 1
  AND no_other_cancer_5y = 1 AND adult = 1 AND prior_obs_365d = 1;

DROP TABLE IF EXISTS #first_platinum;
DROP TABLE IF EXISTS #codesets;
