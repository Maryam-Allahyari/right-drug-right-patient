# Right Drug, Right Patient

**How often are subgroup and biomarker findings in an oncology trial real, and which methods can you trust?**

A simulation study built on a real CDISC reporting pipeline, using Roche's open-source R tools.

**Progress:** planning. The statistical analysis plan is written ([docs/sap.md](docs/sap.md)); The repo will be updater with code later on.

## Why this question
Trials report a forest plot: one treatment effect per subgroup (sex, age, region, biomarker level...). With many subgroups, some will look different by chance alone, and acting on a false subgroup result means giving a drug to patients it doesn't help, or withholding it from patients it does. Biomarker-driven treatment choice (personalised healthcare) depends on telling real effects from noise.

## Methode
In a real trial nobody knows the true effect in each subgroup, so there is no way to check whether a method got it right. This project simulates trials where the truth is known. Each trial is written out as SDTM, converted to ADaM analysis datasets, and analysed with standard and shrinkage subgroup methods. Repeating this thousands of times shows how often each method finds the truth.

- 4 possible scenarios: no effect; same effect for everyone; a real effect in one subgroup; an effect that changes with a continuous biomarker
- 2 trial sizes: 300 and 800 patients
- Tests: standard per-subgroup Cox models, interaction tests, penalised and Bayesian shrinkage, continuous vs dichotomised biomarker
- Measures: false positives, power, accuracy, interval coverage

All data is **simulated**. No real patient data and no real drug have been used; therefore the results and their interpretation cannot be applied to the real-world settings.

## Tools (planned)
R · CDISC SDTM/ADaM · `admiral`, `admiralonco` · `rtables`, `tern` · `bonsaiforest` (brms/Stan) · `survival` · Quarto · Docker · renv

## Results

## Discussion

## Refrences

European Medicines Agency. Guideline on the investigation of subgroups in confirmatory clinical trials. EMA/CHMP/539146/2013. Published 2019. Accessed October 20, 2026. https://www.ema.europa.eu/en/investigation-subgroups-confirmatory-clinical-trials-scientific-guideline
Bornkamp B, Best N, Bratton D, Jelizarow M, Mühlemann N, McMurray J. Bayesian shrinkage methods for routine estimation of subgroup treatment effects [presentation]. European Medicines Agency; 2025. (I couldn’t read the title slide, so check the exact title, event and date there.)
CDISC. Study Data Tabulation Model (SDTM). Accessed October 20, 2026. https://www.cdisc.org/standards/foundational/sdtm
CDISC. Analysis Data Model (ADaM). Accessed October 20, 2026. https://www.cdisc.org/standards/foundational/adam
R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing; 2024. https://www.R-project.org/
Straub B, Bundfuss S, Cascone A, et al. admiral: ADaM in R Asset Library. R package version 1.3.1. 2025. https://CRAN.R-project.org/package=admiral
Bundfuss S, Jain A, Nguyen V, et al. admiralonco: Oncology Extension Package for ADaM in ‘R’ Asset Library. R package version 1.4.0. 2026. https://CRAN.R-project.org/package=admiralonco
Becker G, Waddell A, Garolini D, de la Rua E, Zhu J. rtables: Reporting Tables. R package version 0.6.17. 2026. https://CRAN.R-project.org/package=rtables
Zhu J, Sabanés Bové D, Stoilova J, et al. tern: Create Common TLGs Used in Clinical Trials. R package version 0.9.10. 2025. https://CRAN.R-project.org/package=tern
Vázquez Rabuñal M, Sabanés Bové D, Wolbers M. bonsaiforest: Shrinkage Based Forest Plots. R package version 0.1.1. 2024. https://CRAN.R-project.org/package=bonsaiforest
Therneau T. A Package for Survival Analysis in R. R package version 3.8-6. 2026. https://CRAN.R-project.org/package=survival
Bürkner PC. brms: an R package for Bayesian multilevel models using Stan. J Stat Softw. 2017;80(1):1–28. doi:10.18637/jss.v080.i01
Stan Development Team. RStan: the R interface to Stan. R package version 2.32.7. 2025. https://mc-stan.org/
Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. J Stat Softw. 2010;33(1):1–22. doi:10.18637/jss.v033.i01
Sjoberg D, Baillie M, Fruechtenicht C, Haesendonckx S, Treis T. ggsurvfit: Flexible Time-to-Event Figures. R package version 1.2.0. 2025. https://CRAN.R-project.org/package=ggsurvfit
Gower-Page C, Martin K. diffdf: Dataframe Difference Tool. R package version 1.1.2. 2025. https://CRAN.R-project.org/package=diffdf
Paul-Christian Bürkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1-28. doi:10.18637/jss.v080.i01
Friedman J, Hastie T, Tibshirani R (2010). “Regularization Paths for Generalized Linear Models via Coordinate Descent.” Journal of Statistical Software, 33(1), 1–22. doi:10.18637/jss.v033.i01.
