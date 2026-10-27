# Statistical analysis plan (simulation study)

*Version 2, 2026-10-27. Version 1 was written on 2026-10-19; what changed and why is listed at the end.*

## Question
In a first-line trial in advanced non-small cell lung cancer (NSCLC), comparing a new drug with platinum-based chemotherapy, how often are subgroup and biomarker findings real, and which estimation methods can be trusted?

## Simulated trial
- Patients with advanced (stage IV) NSCLC and no previous treatment for advanced disease
- 1:1 randomisation, **Drug X vs platinum-based chemotherapy** (cisplatin + paclitaxel), N = **300** or **800**
- Accrual uniform over 18 months; data cut-off 30 months after the first patient is randomised
- Tumour scans every 6 weeks (±3 days) from randomisation until progression or cut-off
- Random loss to follow-up: 5% per year
- **Primary endpoint: overall survival (OS)**, time from randomisation to death from any cause
- Secondary endpoints: progression-free survival (PFS), scan-based (progression is found at the first scan after it happens; death before a detected progression counts as a PFS event), and objective response (RECIST-style, investigator)

Baseline covariates (independent of each other):

| Variable | Distribution | Prognostic effect (hazard ratio, HR) |
|---|---|---|
| Sex | 40% female, 60% male | none |
| Age | Normal(63, 9), limited to 30–85; subgroups < 65 / ≥ 65 | 1.2 for ≥ 65 |
| Region | North America 30%, Europe 40%, Asia 30% | none |
| ECOG performance status | 0: 40%, 1: 60% | 1.5 for ECOG 1 |
| Histology | non-squamous 70%, squamous 30% | 1.2 for squamous |
| Biomarker | log-normal, median 5, sdlog 0.6; subgroup "high" = > 5 (fixed assay cut-off) | none |

Death and progression are simulated as two competing exponential times. Both hazards are multiplied by the same prognostic and treatment effects, so OS and PFS have the same true hazard ratio. Survival is followed up after progression (SDTM SS domain), so patients alive at the end are censored at the cut-off or when they were lost to follow-up, not at their last scan.

Control arm calibration: median OS about 11 months and median PFS about 5.5 months in the control arm as a whole, close to the platinum-doublet control arms of recent first-line NSCLC trials (median OS 10.7–11.3 months, PFS 4.8–4.9 months; Paz-Ares et al., 2018; Gadgeel et al., 2020). For the best-prognosis group (age < 65, ECOG 0, non-squamous) the medians are about 16 months (OS) and 7.5 months (PFS). Objective response rate in the control arm is about 25%.

## Scenarios (treatment effect on death and progression)
| | Truth |
|---|---|
| S1 Null | HR 1 for everyone |
| S2 Homogeneous | HR 0.70 for everyone |
| S3 One true subgroup | HR 0.55 if ECOG 0, HR 1.00 if ECOG 1 |
| S4 Continuous biomarker | log HR = log(0.70) − 0.5 · log(biomarker / 5) |

True subgroup HRs (the target for every method) are computed by fitting the same Cox models to one very large simulated trial (200,000 patients) per scenario. These are marginal HRs within each subgroup, so they differ slightly from the HRs above, which are conditional on the other covariates.

## Methods (all on OS)
Forest-plot subgroups: sex, age group, region, ECOG, histology, biomarker high/low (13 subgroups, 6 variables).

- **M1 Naive:** a separate Cox model per subgroup (standard forest plot)
- **M2 Interaction tests:** likelihood-ratio test of treatment × variable, one per variable
- **M3 Ridge, M4 Lasso:** `bonsaiforest::elastic_net`
- **Biomarker:** treatment × biomarker interaction with (B1) the fixed cut-off at 5, (B2) linear in log biomarker, (B3) natural spline of log biomarker (3 df)

## Performance measures
- **False heterogeneity rate** (S1, S2): share of trials with ≥ 1 interaction p < 0.05; share with ≥ 1 naive subgroup "significant in the wrong direction" (95% CI excludes 1 on the side opposite to the truth, or excludes 1 when the true HR is 1)
- **Power** (S3: ECOG interaction p < 0.05; S4: B1, B2, B3 p < 0.05)
- **Accuracy:** root mean squared error (RMSE) of the subgroup log HRs against the true subgroup log HRs (M1, M3, M4)
- **Coverage** of the naive 95% CIs
- Each measure is reported with its Monte Carlo standard error (MCSE), so the precision of the simulation itself is visible

## Replications
- 500 trials per scenario × size (4 × 2 × 500 = 4,000 trials). With 500 trials, the MCSE of a proportion is at most 2.2 percentage points.
- **Pipeline validation:** 20 trials per scenario at N = 300 are written out as SDTM (DM, DS, EX, RS, SC, SS), derived to ADaM (ADSL, ADRS, ADTTE) with `admiral`/`admiralonco`, and the derived OS and PFS times and censoring are compared with the values computed straight from the simulation. If they agree exactly, the large study uses the analysis-level data directly.
- **Showcase trial:** one trial (S4, N = 800) is taken through the full pipeline to clinical study report-style tables and figures with `rtables`/`tern`.

## Link to the real-world part
The control arm of this trial is the comparator for the real-world cohort (`omop/`): patients with NSCLC starting platinum-based chemotherapy in routine care. Both are put in the same ADaM structure (ADSL, ADTTE with `PARAMCD = "OS"`), so the same code gives baseline characteristics (age, sex), Kaplan–Meier curves and OS at 12, 24 and 30 months for both.

## Changes from version 1
Version 1 (2026-10-19) described a generic oncology trial with PFS as the primary endpoint. Version 2 was written before the study code in this repository was run. A first exploratory run of version 1 was done locally; its results are not reported.

| Change | Reason |
|---|---|
| The trial is now first-line advanced NSCLC with platinum-based chemotherapy as control | Matches the real-world cohort, which is NSCLC patients receiving platinum-based chemotherapy |
| OS is the primary endpoint; PFS and response are secondary | OS is the only endpoint that health records can measure, so it is the one that links the two parts |
| Death and progression are simulated as competing times with the same treatment effect | Gives OS a known true HR; in version 1 the drug acted on progression only |
| Survival follow-up after progression added (SS domain) | In version 1, patients alive at the end were censored at their last scan, which stops at progression; that shortens OS and is wrong once OS is the primary endpoint |
| Prior treatment lines replaced by histology (squamous / non-squamous) | Everyone is untreated in a first-line trial |
| Sex 40/60, ECOG 0 40%, control medians set from first-line NSCLC trials | More realistic for this disease |
| Bayesian horseshoe (`bonsaiforest::horseshoe`) moved to a later extension | Long computing time; ridge and lasso already show the trade-off |
| MCSE added to all performance measures | Recommended for simulation studies (Morris et al., 2019) |

## References
- Gadgeel S, Rodríguez-Abreu D, Speranza G, et al. Updated analysis from KEYNOTE-189: pembrolizumab or placebo plus pemetrexed and platinum for previously untreated metastatic nonsquamous non-small-cell lung cancer. *J Clin Oncol*. 2020;38(14):1505–1517. doi:10.1200/JCO.19.03136
- Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Stat Med*. 2019;38(11):2074–2102. doi:10.1002/sim.8086
- Paz-Ares L, Luft A, Vicente D, et al. Pembrolizumab plus chemotherapy for squamous non-small-cell lung cancer. *N Engl J Med*. 2018;379(21):2040–2051. doi:10.1056/NEJMoa1810865
