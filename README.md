# Right Drug, Right Patient

**How far can we trust what a cancer trial tells us?**

This project looks at that question from two sides:

1. **Inside the trial:** how often are subgroup and biomarker findings real, and which methods can you trust? A simulation study built on a CDISC reporting pipeline, using Roche's open-source R tools.
2. **Outside the trial** *(in progress)*: do patients treated in routine care resemble the patients in the trial, and how does their survival compare? A real-world cohort built from synthetic electronic health records in the OMOP Common Data Model.

**Progress:**
- Part 1: the statistical analysis plan is written ([docs/sap.md](docs/sap.md)). Code will be added step by step.
- Part 2: the OMOP database and tools are set up, on branch [`omop-comparator`](../../tree/omop-comparator).

All data are **simulated or synthetic**: there are no real patients and no real drug. The results describe how the methods behave in data like these; they say nothing about the effect of any real treatment.

## Part 1: subgroup findings inside the trial

### Why this question
Trials report a forest plot: one treatment effect per subgroup (sex, age, region, biomarker level…). With many subgroups, some will look different by chance alone, and acting on a false subgroup result means giving a drug to patients it doesn't help, or withholding it from patients it does. Biomarker-driven treatment choice (personalised healthcare) depends on telling real effects from noise. Regulators and methodologists warn against over-reading subgroups (EMA, 2019; Kent et al., 2020), and shrinkage methods that pull extreme subgroup estimates towards the overall effect have been proposed as a more reliable alternative (Wolbers et al., 2024).

### Method
In a real trial nobody knows the true effect in each subgroup, so there is no way to check whether a method got it right. This project simulates trials where the truth is known (Morris et al., 2019). Each trial is written out as SDTM, converted to ADaM analysis datasets, and analysed with standard and shrinkage subgroup methods. Repeating this thousands of times shows how often each method finds the truth, and how often it raises a false alarm or misses a real effect.

- 4 scenarios: no effect; the same effect for everyone; a real effect in one subgroup; an effect that changes with a continuous biomarker
- 2 trial sizes: 300 and 800 patients
- Methods: standard per-subgroup Cox models, interaction tests, penalised and Bayesian shrinkage (`bonsaiforest`), continuous vs dichotomised biomarker
- Measures: false positives, false negatives, accuracy, interval coverage

### Tools
R · CDISC SDTM/ADaM · `admiral`, `admiralonco` · `rtables`, `tern` · `bonsaiforest` (`brms`/Stan) · `survival` · Quarto · Docker · renv

## Part 2: real-world patients outside the trial *(in progress)*

### Why this question
Trial patients are selected: they must be fit enough, meet strict criteria and agree to take part. Patients treated in routine care are often older and sicker, so a trial's results may not carry over to them. Real-world data from health records can show who those patients are and how they fare, but only if the cohort is defined carefully and the data are checked.

### Method
Synthetic patient records are generated with Synthea (Walonoski et al., 2018) and converted to the OMOP Common Data Model, the standard used by the OHDSI community (OHDSI, 2021). A cohort of patients starting cancer treatment is defined to mirror the trial's control arm, with a clear start of follow-up and censoring rules, following the ideas of target trial emulation (Hernán & Robins, 2016). The cohort is built twice, in hand-written SQL and with the OHDSI cohort tools, and the two versions are compared as a quality check. Overall survival is then put in the same ADaM format as the trial, so both can be analysed with the same code.

The two data sources were generated independently and never designed to match, so any difference between them says nothing about real patients. The aim is a correct, reproducible and quality-checked pipeline, and a clear account of what would make such a comparison valid with real data.

### Tools (planned)
Synthea · OMOP CDM v5.4 · PostgreSQL / SQL · OHDSI tools (Achilles, DataQualityDashboard, ATLAS, CirceR, CohortGenerator, FeatureExtraction) · SNOMED CT, RxNorm · R

## References

**Papers and guidelines**
- European Medicines Agency. Guideline on the investigation of subgroups in confirmatory clinical trials. EMA/CHMP/539146/2013. 2019. https://www.ema.europa.eu/en/investigation-subgroups-confirmatory-clinical-trials-scientific-guideline
- Hernán MA, Robins JM. Using big data to emulate a target trial when a randomized trial is not available. *Am J Epidemiol*. 2016;183(8):758–764. doi:10.1093/aje/kwv254
- Kent DM, Paulus JK, van Klaveren D, et al. The Predictive Approaches to Treatment effect Heterogeneity (PATH) Statement. *Ann Intern Med*. 2020;172(1):35–45. doi:10.7326/M18-3667
- Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Stat Med*. 2019;38(11):2074–2102. doi:10.1002/sim.8086
- Observational Health Data Sciences and Informatics. *The Book of OHDSI*. 2021. https://ohdsi.github.io/TheBookOfOhdsi/
- Walonoski J, Kramer M, Nichols J, et al. Synthea: an approach, method, and software mechanism for generating synthetic patients and the synthetic electronic health care record. *J Am Med Inform Assoc*. 2018;25(3):230–238. doi:10.1093/jamia/ocx079
- Wolbers M, Vázquez Rabuñal M, Li K, Rufibach K, Sabanés Bové D. Using shrinkage methods to estimate treatment effects in overlapping subgroups in randomized clinical trials with a time-to-event endpoint. *arXiv* [preprint]. 2024. doi:10.48550/arXiv.2407.11729
- Wolbers M, Pedrera Gómez M, Ocampo A, Gravestock I. Unified implementation and comparison of Bayesian shrinkage methods for treatment effect estimation in subgroups. *arXiv* [preprint]. 2026. doi:10.48550/arXiv.2603.21967

**Standards**
- CDISC. Study Data Tabulation Model (SDTM). https://www.cdisc.org/standards/foundational/sdtm
- CDISC. Analysis Data Model (ADaM). https://www.cdisc.org/standards/foundational/adam

**Software**
- Becker G, Waddell A, Garolini D, de la Rua E, Zhu J. *rtables: Reporting Tables*. R package version 0.6.17. 2026. https://CRAN.R-project.org/package=rtables
- Bundfuss S, Jain A, Nguyen V, et al. *admiralonco: Oncology Extension Package for ADaM in 'R' Asset Library*. R package version 1.4.0. 2026. https://CRAN.R-project.org/package=admiralonco
- Bürkner PC. brms: an R package for Bayesian multilevel models using Stan. *J Stat Softw*. 2017;80(1):1–28. doi:10.18637/jss.v080.i01
- Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33(1):1–22. doi:10.18637/jss.v033.i01
- Gower-Page C, Martin K. *diffdf: Dataframe Difference Tool*. R package version 1.1.2. 2025. https://CRAN.R-project.org/package=diffdf
- R Core Team. *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing; 2024. https://www.R-project.org/
- Sjoberg D, Baillie M, Fruechtenicht C, Haesendonckx S, Treis T. *ggsurvfit: Flexible Time-to-Event Figures*. R package version 1.2.0. 2025. https://CRAN.R-project.org/package=ggsurvfit
- Stan Development Team. *RStan: the R interface to Stan*. R package version 2.32.7. 2025. https://mc-stan.org/
- Straub B, Bundfuss S, Cascone A, et al. *admiral: ADaM in R Asset Library*. R package version 1.3.1. 2025. https://CRAN.R-project.org/package=admiral
- Therneau T. *A Package for Survival Analysis in R*. R package version 3.8-6. 2026. https://CRAN.R-project.org/package=survival
- Vázquez Rabuñal M, Sabanés Bové D, Wolbers M. *bonsaiforest: Shrinkage Based Forest Plots*. R package version 0.1.1. 2024. https://CRAN.R-project.org/package=bonsaiforest
- Zhu J, Sabanés Bové D, Stoilova J, et al. *tern: Create Common TLGs Used in Clinical Trials*. R package version 0.9.10. 2025. https://CRAN.R-project.org/package=tern
