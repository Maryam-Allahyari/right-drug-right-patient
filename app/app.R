# Interactive version of the showcase trial (S4, 800 patients): the same
# ADaM datasets as the CSR-style tables in trial/results/showcase, explored
# with teal. Pick the endpoint, strata, subgroups and time points yourself
# instead of reading fixed tables. Every output can show the R code that
# made it ("Show R code"), so results stay reproducible.
#
# Before running:
#   - trial/scripts/05_showcase.R run (it saves adsl.rds, adtte.rds, adrs.rds)
#   - install.packages("teal.modules.clinical")   # brings teal and tern
#
# From the repo root:
#   Rscript -e 'shiny::runApp("app")'

library(teal.modules.clinical)

# shiny runs the app from its own folder, so the data is one level up
showcase <- file.path("..", "trial", "results", "showcase")

# CDISC-style labels for the variables the modules show. teal uses the
# labels to build its selectors and table rows, and fails without them;
# columns not listed here are labelled with their own name
labels <- c(
  STUDYID = "Study Identifier", USUBJID = "Unique Subject Identifier",
  ARM = "Description of Planned Arm", ARMCD = "Planned Arm Code",
  AGE = "Age", AGEGR1 = "Age Group", SEX = "Sex", REGION1 = "Geographic Region",
  ECOG = "ECOG Performance Status", HISTOL = "Histology",
  BMRK = "Biomarker Value", BMRKGR1 = "Biomarker Group (cut-off 5)",
  PARAMCD = "Parameter Code", PARAM = "Parameter", AVAL = "Analysis Value",
  AVALU = "Analysis Value Unit", AVALC = "Analysis Value (C)", CNSR = "Censor",
  EVNTDESC = "Event or Censoring Description", AVISIT = "Analysis Visit"
)

data <- within(teal_data(), {
  library(dplyr)
  add_labels <- function(d) {
    lab <- setNames(names(d), names(d))
    known <- intersect(names(d), names(labels))
    lab[known] <- labels[known]
    formatters::var_labels(d) <- lab
    d
  }
  ADSL <- add_labels(readRDS(file.path(showcase, "adsl.rds")))
  # time in months reads better in KM plots and tables than days
  ADTTE <- readRDS(file.path(showcase, "adtte.rds")) %>%
    mutate(AVAL = AVAL / 30.4375, AVALU = "MONTHS") %>%
    add_labels()
  ADRS <- add_labels(readRDS(file.path(showcase, "adrs.rds")))
}, showcase = showcase, labels = labels)
join_keys(data) <- default_cdisc_join_keys[c("ADSL", "ADTTE", "ADRS")]

arm <- choices_selected(c("ARM", "ARMCD"), "ARM")
arm_ref_comp <- list(ARM = list(ref = "Platinum chemotherapy", comp = "Drug X"),
                     ARMCD = list(ref = "CHEMO", comp = "DRUGX"))
baseline <- c("SEX", "AGEGR1", "REGION1", "ECOG", "HISTOL", "BMRKGR1")
tte_param <- choices_selected(value_choices(data[["ADTTE"]], "PARAMCD", "PARAM"), "OS")

app <- init(
  data = data,
  modules = modules(
    tm_t_summary(
      label = "Demographics", dataname = "ADSL", arm_var = arm,
      summarize_vars = choices_selected(c("AGE", baseline, "BMRK"), c("AGE", baseline))
    ),
    tm_g_km(
      label = "Kaplan-Meier plot", dataname = "ADTTE", arm_var = arm,
      arm_ref_comp = arm_ref_comp, paramcd = tte_param,
      strata_var = choices_selected(baseline, NULL),
      facet_var = choices_selected(baseline, NULL)
    ),
    tm_t_tte(
      label = "Time-to-event table", dataname = "ADTTE", arm_var = arm,
      arm_ref_comp = arm_ref_comp, paramcd = tte_param,
      strata_var = choices_selected(baseline, NULL),
      time_points = choices_selected(c(6, 12, 24), 12),
      event_desc_var = choices_selected("EVNTDESC", "EVNTDESC", fixed = TRUE)
    ),
    tm_g_forest_tte(
      label = "Forest plot", dataname = "ADTTE", arm_var = arm,
      arm_ref_comp = arm_ref_comp, paramcd = tte_param,
      subgroup_var = choices_selected(baseline, baseline),
      strata_var = choices_selected(baseline, NULL)
    ),
    tm_t_binary_outcome(
      label = "Response", dataname = "ADRS", arm_var = arm,
      arm_ref_comp = arm_ref_comp,
      paramcd = choices_selected(value_choices(data[["ADRS"]], "PARAMCD", "PARAM",
                                               subset = c("BOR", "RSP")), "BOR"),
      strata_var = choices_selected(baseline, NULL)
    )
  )
) |>
  modify_title("NSCLC trial (simulated): Drug X vs platinum chemotherapy") |>
  modify_header("Showcase trial, scenario S4, 800 patients. All data are simulated.")

shinyApp(app$ui, app$server)
