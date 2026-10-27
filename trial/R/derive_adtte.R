# ADTTE: overall survival (OS, primary) and progression-free survival (PFS).
# Start = randomisation. Event and censoring rules are the admiralonco ones:
#   OS : event = death; censored at the last date known alive
#   PFS: event = first PD or death; censored at the last tumour assessment
# Patients with no data at all are censored at randomisation (rand_censor).

library(admiral)
library(admiralonco)
library(dplyr)

derive_adtte <- function(adsl, adrs) {
  src <- list(adsl = adsl, adrs = adrs)
  keys <- get_admiral_option("subject_keys")

  derive_param_tte(
    dataset_adsl = adsl, start_date = RANDDT,
    event_conditions = list(death_event),
    censor_conditions = list(lastalive_censor, rand_censor),
    source_datasets = src,
    set_values_to = exprs(PARAMCD = "OS", PARAM = "Overall Survival")
  ) %>%
    derive_param_tte(
      dataset_adsl = adsl, start_date = RANDDT,
      event_conditions = list(pd_event, death_event),
      censor_conditions = list(lasta_censor, rand_censor),
      source_datasets = src,
      set_values_to = exprs(PARAMCD = "PFS", PARAM = "Progression Free Survival")
    ) %>%
    # AVAL in days, counting the start day as day 1
    derive_vars_duration(new_var = AVAL, start_date = STARTDT, end_date = ADT) %>%
    derive_var_obs_number(by_vars = keys, order = exprs(PARAMCD), check_type = "error") %>%
    # treatment and the subgroup variables, so the TTE data can be analysed on its own
    derive_vars_merged(dataset_add = adsl, by_vars = keys,
                       new_vars = exprs(ARMCD, ARM, TRT01A, AGE, AGEGR1, SEX, REGION1,
                                        ECOG, HISTOL, BMRK, BMRKGR1))
}
