# Right Drug, Right Patient

**How far can we trust what a cancer trial tells us?**

This project looks at that question from two sides, using advanced non-small cell lung cancer (NSCLC) as the example:

1. **Inside the trial:** in a first-line NSCLC trial comparing a new drug with platinum chemotherapy, which subgroup findings can you trust? A simulation study built on a CDISC reporting pipeline, using Roche's open-source R tools.
2. **Outside the trial** *(in progress)*: the patients who get platinum chemotherapy in routine care. Do they look like that trial's control arm, and how long do they live? A real-world cohort built from synthetic health records in the OMOP Common Data Model.

**Progress:**
- Part 1: done. The plan ([docs/sap.md](docs/sap.md)), code ([trial/](trial/)) and results are below.
- Part 2: synthetic patients generated and profiled, mapping to OMOP documented, data being loaded into the OMOP database ([omop/](omop/)). The cohort comes next.

All data are **simulated or synthetic**: there are no real patients and no real drug. The results describe how the methods behave in data like these; they say nothing about the effect of any real treatment.

## How the two parts fit together

| | Part 1: simulated trial | Part 2: real-world cohort |
|---|---|---|
| Patients | first-line advanced NSCLC | NSCLC patients in synthetic health records |
| Treatment compared | Drug X vs **platinum chemotherapy** (cisplatin + paclitaxel) | **platinum chemotherapy** (cisplatin, paclitaxel, etoposide) |
| Main endpoint | overall survival (OS) | OS (health records have no tumour scans, so no progression) |
| Data standard | CDISC SDTM → ADaM | OMOP CDM → the same ADaM structure (ADSL, ADTTE) |

The trial's control arm and the real-world cohort get the same treatment and end up in the same ADaM format, so one piece of code can describe both: who the patients are (age, sex) and how long they live (Kaplan–Meier curves, OS at 12, 24 and 30 months).

## Part 1: subgroup findings inside the trial

### Why this question
Trials report a forest plot: one treatment effect per subgroup (sex, age, region, biomarker level…). With many subgroups, some will look different by chance alone. Acting on a false subgroup result means giving a drug to patients it doesn't help, or withholding it from patients it does, and biomarker-driven treatment (personalised healthcare) depends on telling real effects from noise. Regulators and methodologists warn against over-reading subgroups (EMA, 2019; Kent et al., 2020). Shrinkage methods, which pull extreme subgroup estimates towards the overall effect, have been proposed as a more reliable alternative (Wolbers et al., 2024).

### Design
In a real trial nobody knows the true effect in each subgroup, so there's no way to check whether a method got it right. Here the trials are simulated, so the truth is known (Morris et al., 2019). The plan was written before the code was run ([docs/sap.md](docs/sap.md), version 2; the changes from version 1 are listed at its end).

- **Trial:** first-line advanced NSCLC, 1:1 Drug X vs cisplatin + paclitaxel, 300 or 800 patients, 30 months of follow-up. The control arm is set up like published platinum-chemotherapy arms (median OS about 11 months; Paz-Ares et al., 2018; Gadgeel et al., 2020).
- **Primary endpoint:** OS. Progression-free survival (PFS) and tumour response are secondary.
- **Subgroups:** sex, age (< 65 / ≥ 65), region, ECOG performance status (Eastern Cooperative Oncology Group), histology (squamous or not) and a biomarker with an assay cut-off: 6 variables, 13 subgroups.
- **Scenarios:** S1 the drug does nothing; S2 the same effect for everyone (hazard ratio, HR 0.70); S3 an effect only in ECOG 0 patients (HR 0.55 vs 1.00); S4 an effect that grows with the biomarker.
- **Methods:** a separate Cox model per subgroup (the standard forest plot), treatment-by-subgroup interaction tests, ridge and lasso shrinkage (`bonsaiforest`), and the biomarker modelled three ways (cut-off, linear on the log scale, spline).
- **Size:** 500 trials per scenario and trial size (4,000 trials). Every result has a Monte Carlo standard error (MCSE) of at most 2.2 percentage points.

### Checking the CDISC pipeline
80 simulated trials (24,000 patients) were written out as SDTM (DM, DS, EX, RS, SC, SS), derived to ADaM (ADSL, ADRS, ADTTE) with `admiral` and `admiralonco`, and compared patient by patient with the OS and PFS known from the simulation. **There were no differences** in survival times or censoring, so the large study works on the analysis data directly ([pipeline_vs_truth.csv](trial/results/pipeline_vs_truth.csv)).

### Key findings
- **When the drug does nothing, a third of trials (33%) still show at least one "significant" subgroup** in the standard forest plot, and 27% show a "significant" interaction. With 800 patients instead of 300 it's the same (32%): the problem is the number of looks, not the trial size.
- **When the drug works equally well for everyone, interaction tests still suggest a difference in 25–29% of trials.**
- **A real subgroup effect is easy to miss and easy to dilute.** When only ECOG 0 patients benefit, the overall result (HR 0.81) is significant in only 30% of trials with 300 patients. The ECOG interaction is found in 47% of trials with 300 patients and 87% with 800, and in about 1 in 5 trials one of the four variables with no real effect is flagged as well.
- **Shrinkage makes subgroup estimates about a third more accurate when there are no real differences** (root mean squared error, RMSE, 0.21 → 0.14 with 300 patients).
- **But ridge shrinkage erases a real difference.** In S3 it gives ECOG 0 and ECOG 1 almost the same HR (0.83 and 0.83 with 800 patients; truth 0.55 and 1.01), and more patients don't fix it: with 800 patients ridge is *less* accurate than the standard forest plot. Lasso keeps much of the difference (0.64 vs 0.94) and is the most accurate method in S3 and S4 with 800 patients.
- **Keeping the biomarker continuous finds a real biomarker effect far more often** than the high/low cut-off (51% vs 33% of trials with 300 patients, 93% vs 79% with 800), with the same false-positive rate (4–7%).
- The standard per-subgroup confidence intervals have correct coverage (95–96%). The estimates themselves are fine; the problem is reading 13 of them at once.

![S3: shrinkage and a real subgroup effect](trial/results/fig_s3_signal.png)

All measures with their MCSEs: [trial/results/summary.csv](trial/results/summary.csv). Figures: [accuracy](trial/results/fig_rmse.png), [biomarker models](trial/results/fig_biomarker_power.png).

### One trial, reported like a real study
One simulated trial (S4, 800 patients) goes through the full pipeline to clinical study report-style outputs with `rtables` and `tern` ([trial/results/showcase/](trial/results/showcase/)):

| Output | Result |
|---|---|
| [Table 14.1.1](trial/results/showcase/t_14_1_1_demographics.txt) Demographics | 800 patients, median age 64, 62% male |
| [Table 14.2.1](trial/results/showcase/t_14_2_1_os.txt), [Figure 14.2.1](trial/results/showcase/f_14_2_1_km_os.png) OS | median 15.4 vs 11.1 months, HR 0.72 (95% CI 0.60–0.86) |
| [Table 14.2.2](trial/results/showcase/t_14_2_2_pfs.txt) PFS | median 7.2 vs 5.2 months, HR 0.66 |
| [Table 14.2.3](trial/results/showcase/t_14_2_3_response.txt) Response | 37.0% vs 24.5% |
| [Figure 14.2.2](trial/results/showcase/f_14_2_2_forest_os.png) Forest plot | women 0.82, men 0.65, although the true effect doesn't depend on sex |
| [Figure 14.2.3](trial/results/showcase/f_14_2_3_hr_by_biomarker.png) HR across the biomarker | the estimated curve follows the true one; interaction p = 0.02 |

### Running it
From the repo root, in order (R 4.4; packages listed at the top of each script):

```
Rscript trial/scripts/01_validate_pipeline.R   # SDTM -> ADaM vs the simulated truth (~3 min)
Rscript trial/scripts/02_truth.R               # true subgroup HRs
Rscript trial/scripts/03_study.R               # 4,000 trials, in parallel (~8 min on 16 cores)
Rscript trial/scripts/04_summarise.R           # performance measures, MCSEs, figures
Rscript trial/scripts/05_showcase.R            # one trial to tables and figures
```

## Part 2: real-world patients outside the trial *(in progress)*

### Why this question
Trial patients are selected: they must be fit enough, meet strict criteria and agree to take part. Patients treated in routine care are often older and sicker, so a trial's results may not carry over to them. Health records can show who those patients are and how they fare, but only if the cohort is defined carefully and the data are checked.

### Method
Synthetic patient records are generated with Synthea (Walonoski et al., 2018) and converted to the OMOP Common Data Model, the standard used by the OHDSI community (OHDSI, 2021). The cohort is NSCLC patients starting platinum chemotherapy, defined to mirror the trial's control arm, with a clear start of follow-up and censoring rules, following the ideas of target trial emulation (Hernán & Robins, 2016). It is built twice, in hand-written SQL and with the OHDSI cohort tools, and the two versions are compared as a quality check. OS is then put in the same ADaM format as the trial.

The trial and the synthetic records were generated independently and never designed to match. Synthea's lung cancer patients, for example, live much longer than real ones (see [omop/README.md](omop/README.md)). So any difference between the two says nothing about real patients. The aim is a correct, reproducible and quality-checked pipeline, and a clear account of what would make such a comparison valid with real data.

Done so far: choice of cancer from two pilot populations, the full population (835 lung cancer patients), a WhiteRabbit profile of the source data, the Synthea → OMOP mapping in Rabbit-in-a-Hat, and the OMOP vocabularies. Details in [omop/README.md](omop/README.md).

## Repository layout

```
docs/sap.md          analysis plan for Part 1
trial/R/             simulator, SDTM writer, ADaM derivations, analysis methods
trial/scripts/       01-05, run in order
trial/results/       results, figures, showcase tables and SDTM
omop/                Part 2: Synthea, WhiteRabbit, mapping, OMOP load
```

## Tools
**Part 1:** R · CDISC SDTM/ADaM · `admiral`, `admiralonco` · `rtables`, `tern` · `bonsaiforest` · `survival`, `ggsurvfit`, `ggplot2`

**Part 2:** Synthea · OMOP CDM v5.4 · PostgreSQL / SQL · OHDSI tools (WhiteRabbit, Rabbit-in-a-Hat, ETL-Synthea; planned: Achilles, DataQualityDashboard, ATLAS, CirceR, CohortGenerator) · SNOMED CT, RxNorm · R

## References

**Papers and guidelines**
- European Medicines Agency. Guideline on the investigation of subgroups in confirmatory clinical trials. EMA/CHMP/539146/2013. 2019. https://www.ema.europa.eu/en/investigation-subgroups-confirmatory-clinical-trials-scientific-guideline
- Gadgeel S, Rodríguez-Abreu D, Speranza G, et al. Updated analysis from KEYNOTE-189: pembrolizumab or placebo plus pemetrexed and platinum for previously untreated metastatic nonsquamous non-small-cell lung cancer. *J Clin Oncol*. 2020;38(14):1505–1517. doi:10.1200/JCO.19.03136
- Hernán MA, Robins JM. Using big data to emulate a target trial when a randomized trial is not available. *Am J Epidemiol*. 2016;183(8):758–764. doi:10.1093/aje/kwv254
- Kent DM, Paulus JK, van Klaveren D, et al. The Predictive Approaches to Treatment effect Heterogeneity (PATH) Statement. *Ann Intern Med*. 2020;172(1):35–45. doi:10.7326/M18-3667
- Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Stat Med*. 2019;38(11):2074–2102. doi:10.1002/sim.8086
- Observational Health Data Sciences and Informatics. *The Book of OHDSI*. 2021. https://ohdsi.github.io/TheBookOfOhdsi/
- Paz-Ares L, Luft A, Vicente D, et al. Pembrolizumab plus chemotherapy for squamous non-small-cell lung cancer. *N Engl J Med*. 2018;379(21):2040–2051. doi:10.1056/NEJMoa1810865
- Walonoski J, Kramer M, Nichols J, et al. Synthea: an approach, method, and software mechanism for generating synthetic patients and the synthetic electronic health care record. *J Am Med Inform Assoc*. 2018;25(3):230–238. doi:10.1093/jamia/ocx079
- Wolbers M, Vázquez Rabuñal M, Li K, Rufibach K, Sabanés Bové D. Using shrinkage methods to estimate treatment effects in overlapping subgroups in randomized clinical trials with a time-to-event endpoint. *arXiv* [preprint]. 2024. doi:10.48550/arXiv.2407.11729
- Wolbers M, Pedrera Gómez M, Ocampo A, Gravestock I. Unified implementation and comparison of Bayesian shrinkage methods for treatment effect estimation in subgroups. *arXiv* [preprint]. 2026. doi:10.48550/arXiv.2603.21967

**Standards**
- CDISC. Study Data Tabulation Model (SDTM). https://www.cdisc.org/standards/foundational/sdtm
- CDISC. Analysis Data Model (ADaM). https://www.cdisc.org/standards/foundational/adam
- OHDSI. OMOP Common Data Model v5.4. https://ohdsi.github.io/CommonDataModel/

**Software**
- Becker G, Waddell A, Garolini D, de la Rua E, Zhu J. *rtables: Reporting Tables*. R package version 0.6.17. 2026. https://CRAN.R-project.org/package=rtables
- Bundfuss S, Jain A, Nguyen V, et al. *admiralonco: Oncology Extension Package for ADaM in 'R' Asset Library*. R package version 1.4.0. 2026. https://CRAN.R-project.org/package=admiralonco
- Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33(1):1–22. doi:10.18637/jss.v033.i01
- R Core Team. *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing; 2024. https://www.R-project.org/
- Sjoberg D, Baillie M, Fruechtenicht C, Haesendonckx S, Treis T. *ggsurvfit: Flexible Time-to-Event Figures*. R package version 1.2.0. 2025. https://CRAN.R-project.org/package=ggsurvfit
- Straub B, Bundfuss S, Cascone A, et al. *admiral: ADaM in R Asset Library*. R package version 1.3.1. 2025. https://CRAN.R-project.org/package=admiral
- Therneau T. *A Package for Survival Analysis in R*. R package version 3.8-6. 2026. https://CRAN.R-project.org/package=survival
- Vázquez Rabuñal M, Sabanés Bové D, Wolbers M. *bonsaiforest: Shrinkage Based Forest Plots*. R package version 0.1.1. 2024. https://CRAN.R-project.org/package=bonsaiforest
- Zhu J, Sabanés Bové D, Stoilova J, et al. *tern: Create Common TLGs Used in Clinical Trials*. R package version 0.9.10. 2025. https://CRAN.R-project.org/package=tern
