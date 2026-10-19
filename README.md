# Right Drug, Right Patient

**How often are subgroup and biomarker findings in an oncology trial real, and which methods can you trust?**

A simulation study built on a real CDISC reporting pipeline, using Roche's open-source R tools.

**Status:** planning. The statistical analysis plan is written ([docs/sap.md](docs/sap.md)); code comes next.

## Why this question
Trials report a forest plot: one treatment effect per subgroup (sex, age, region, biomarker level...). With many subgroups, some will look different by chance alone, and acting on a false subgroup result means giving a drug to patients it doesn't help, or withholding it from patients it does. Biomarker-driven treatment choice (personalised healthcare) depends on telling real effects from noise.

## Approach
Because the truth is never known in a real trial, this project simulates trials where it is:

```
simulate trials (known truth) -> SDTM (DM, DS, EX, RS) -> ADaM (ADSL, ADRS, ADTTE)
  -> standard and shrinkage subgroup methods -> repeat many times -> compare with the truth
```

- 4 scenarios: no effect; same effect for everyone; a real effect in one subgroup; an effect that changes with a continuous biomarker
- 2 trial sizes: 300 and 800 patients
- Methods: standard per-subgroup Cox models, interaction tests, penalised and Bayesian shrinkage, continuous vs dichotomised biomarker
- Measures: false positives, power, accuracy, interval coverage

All data is **simulated**. No real patient data and no real drug.

## Tools (planned)
R · CDISC SDTM/ADaM · `admiral`, `admiralonco` · `rtables`, `tern` · `bonsaiforest` (brms/Stan) · `survival` · Quarto · Docker · renv

## Repository layout (planned)
```
R/          functions (simulation, SDTM writer, ADaM derivation, methods)
scripts/    numbered analysis steps
docs/       analysis plan, specs, QC log
outputs/    generated tables and figures (not tracked)
```
