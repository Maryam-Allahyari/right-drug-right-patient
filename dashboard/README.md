# Tableau dashboard

The main results as two Tableau dashboards. `export_tableau.R` writes everything Tableau needs into `data/` as tidy CSV files (one row per number), so the dashboards hold no calculations of their own: if a result changes, rerun the export and refresh the data source in Tableau.

```
Rscript dashboard/export_tableau.R
```

## The data

| File | One row per | Used for |
|---|---|---|
| `study_measures.csv` | scenario × trial size × measure (× method), with MCSE | error rates, power, accuracy (RMSE) of the methods |
| `s3_ecog_signal.csv` | trial size × method × ECOG group | does the method keep the real ECOG difference? |
| `delayed_effect.csv` | scenario × trial size | log-rank vs RMST, average HR |
| `forest_showcase.csv` | subgroup | forest plot of the showcase trial |
| `km_curves.csv` | panel × group × event time | Kaplan–Meier curves (showcase arms; trial control arm vs routine care) |
| `baseline_trial_vs_routine_care.csv` | source × measure | age and sex, trial vs routine care |
| `cohort_attrition.csv` | inclusion rule | the cohort's patient flow |
| `hybrid_control.csv` | borrowing analysis | HR with and without external controls |
| `data_quality_overview.csv` | check category | DataQualityDashboard results |

`study_measures.csv` has a `trials` column: 500 for the main study, 200 for the horseshoe, which was run on a subset. Filter or colour by it so the two aren't compared as if they were the same trials.

## Dashboard 1: Which subgroup findings can you trust?

1. **False alarms** (`study_measures`): `measure` on Filters (keep *At least one misleading subgroup* and *Any interaction test significant*), `scenario_label` on Columns (S1 and S2 only), `value` on Rows as a bar, `n` on Color. Add a reference line at 0.05: the error rate one test is supposed to have.
2. **Accuracy of the methods** (`study_measures`): filter `measure` = *RMSE of subgroup log HR*, `method` on Columns, `value` on Rows as a bar, `scenario_label` and `n` as small multiples. Error bars: calculated fields `value - 1.96 * mcse` and `value + 1.96 * mcse`, drawn as a Gantt bar or as reference bands.
3. **Forest plot of one trial** (`forest_showcase`): `variable` and `level` on Rows, `hr` on Columns with a logarithmic axis, Gantt bar from `lower` with size `upper - lower`, a circle at `hr`, reference line at 1. The tooltip shows patients and events.
4. **Biomarker** (`study_measures`): filter the three *Biomarker effect found* measures, scenario S4, `n` on Color.

Add a filter on trial size (`n`) applied to all sheets.

## Dashboard 2: Do routine-care patients look like the trial's control arm?

1. **Survival** (`km_curves`): `months` on Columns, `survival` on Rows, `group` on Color, mark type Line with the line type set to *Step*. `panel` as a filter (single value) to switch between the showcase trial's arms and the trial control arm vs routine care. `lower` and `upper` can be added as a band with a dual axis.
2. **Who they are** (`baseline_trial_vs_routine_care`): `measure` on Rows, `value` on Columns as a bar, `source` on Color.
3. **Patient flow** (`cohort_attrition`): `step` on Rows in file order, `patients` on Columns as a bar, labels showing `removed`.
4. **Borrowing external controls** (`hybrid_control`): `analysis` on Rows, Gantt bar from `lower` to `upper`, circle at `hr`, logarithmic axis, reference line at 1. Full borrowing from routine care sits far to the right of all the others.

## Publishing

Tableau Public (free) needs an account and Tableau Public Desktop. Build the sheets and dashboards as above, then *File → Save to Tableau Public*. Put the link in the main README next to the results it shows. The published workbook keeps a copy of the data, so republish after rerunning the export.
