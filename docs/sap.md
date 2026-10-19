# Statistical analysis plan (simulation study)

*Written 2026-10-19, before any simulation was run. Changes after this date are listed at the end.*

## Question
In a randomised oncology trial of realistic size, how often are subgroup and biomarker findings real, and which estimation methods can be trusted?

## Simulated trial
- 1:1 randomisation, Drug vs Control, N = **300** or **800**
- Accrual uniform over 18 months; data cut-off 30 months after first patient in
- Tumour scans every 6 weeks (±3 days) from randomisation until progression or cut-off
- Random loss to follow-up: 5% per year
- **Primary endpoint: PFS** (investigator, scan-based: progression is detected at the first scan after it happens). Death before a detected progression is a PFS event. OS is secondary.

Baseline covariates (independent):

| Variable | Distribution | Prognostic effect on progression (HR) |
|---|---|---|
| Sex | 50% F / 50% M | none |
| Age | Normal(62, 10); subgroup < 65 / ≥ 65 | 1.2 for ≥ 65 |
| Region | North America 30%, Europe 40%, Asia 30% | none |
| ECOG | 0: 55%, 1: 45% | 1.6 for ECOG 1 |
| Prior lines | 1: 60%, ≥ 2: 40% | 1.3 for ≥ 2 |
| Biomarker | log-normal, median 5, sdlog 0.6; subgroup "high" = > 5 (fixed assay cut-off) | none |

Control arm median time to progression ≈ 7 months; death before progression at 1% per month; survival after progression exponential, median 10 months.

## Scenarios (treatment effect on progression)
| | Truth |
|---|---|
| S1 Null | HR 1 for everyone |
| S2 Homogeneous | HR 0.70 for everyone |
| S3 One true subgroup | HR 0.55 if ECOG 0, HR 1.00 if ECOG 1 |
| S4 Continuous biomarker | log HR = log(0.70) − 0.5 · log(biomarker / 5) |

True subgroup HRs (the target for every method) are computed by fitting the same Cox models to one very large simulated trial (200,000 patients) per scenario.

## Methods
Forest-plot subgroups: sex, age group, region, ECOG, prior lines, biomarker high/low (13 subgroups, 6 variables).

- **M1 Naive:** a separate Cox model per subgroup (standard forest plot)
- **M2 Interaction tests:** likelihood-ratio test of treatment × variable, one per variable
- **M3 Ridge, M4 Lasso:** `bonsaiforest::elastic_net`
- **M5 Bayesian horseshoe:** `bonsaiforest::horseshoe`, on a subset only (computing time)
- **Biomarker:** treatment × biomarker interaction with (B1) the fixed cut-off at 5, (B2) linear in log biomarker, (B3) natural spline of log biomarker (3 df)

## Performance measures
- **False heterogeneity rate** (S1, S2): share of trials with ≥ 1 interaction p < 0.05; share with ≥ 1 naive subgroup "significant in the wrong direction" (CI excludes 1 on the opposite side of the truth)
- **Power** (S3: ECOG interaction p < 0.05; S4: B1, B2, B3 p < 0.05)
- **Accuracy:** RMSE of the subgroup log HRs vs the true subgroup log HRs (M1, M3, M4, M5)
- **Coverage** of 95% intervals (M1 CIs, M5 credible intervals)

## Replications
- 500 trials per scenario × size for M1–M4 and B1–B3 (4 × 2 × 500 = 4,000 trials)
- M5: 25 trials per scenario at N = 300 (100 fits)
- **Pipeline validation:** 20 trials per scenario at N = 300 are written out as SDTM (DM, DS, EX, RS), derived to ADaM with admiral/admiralonco, and the derived PFS/OS and best response are compared with the values computed straight from the simulation. If they agree exactly, the large study uses the analysis-level data directly.

## Changes after this date
(none yet)
