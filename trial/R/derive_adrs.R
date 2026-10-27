# ADRS: tumour response (RECIST 1.1, investigator).
# Parameters: OVR (per visit), PD, RSP, BOR, CRSP, CBOR, DEATH, LSTA.
# Based on the admiralonco ADRS template. MDIS and the clinical benefit
# parameters are left out, they aren't used downstream.

library(admiral)
library(admiralonco)
library(dplyr)

# ranking used to pick one response if a visit has several records
worst_resp <- function(x) {
  case_when(x == "NE" ~ 1, x == "CR" ~ 2, x == "PR" ~ 3, x == "SD" ~ 4,
            x == "NON-CR/NON-PD" ~ 5, x == "PD" ~ 6, TRUE ~ 0)
}

recist_labels <- exprs(PARCAT1 = "Tumor Response", PARCAT2 = "Investigator",
                       PARCAT3 = "RECIST 1.1")

derive_adrs <- function(rs, adsl) {
  keys <- get_admiral_option("subject_keys")
  adsl_vars <- exprs(RANDDT)

  adrs <- convert_blanks_to_na(rs) %>%
    derive_vars_merged(dataset_add = adsl, new_vars = adsl_vars, by_vars = keys) %>%
    filter(RSEVAL == "INVESTIGATOR", RSTESTCD == "OVRLRESP") %>%
    mutate(PARAMCD = "OVR", PARAM = "Overall Response by Investigator",
           !!!recist_labels) %>%
    derive_vars_dt(dtc = RSDTC, new_vars_prefix = "A",
                   highest_imputation = "D", date_imputation = "last") %>%
    mutate(AVISIT = VISIT, AVALC = RSSTRESC, AVAL = aval_resp(AVALC))

  # ANL01FL: one record per date, after randomisation
  # ANL02FL: records up to and including the first PD
  adrs <- adrs %>%
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = c(keys, exprs(ADT)),
                    order = exprs(worst_resp(AVALC), RSSEQ),
                    new_var = ANL01FL, mode = "last"),
      filter = !is.na(AVAL) & ADT >= RANDDT
    ) %>%
    derive_var_relative_flag(
      by_vars = keys, order = exprs(ADT, RSSEQ), new_var = ANL02FL,
      condition = AVALC == "PD", mode = "first", selection = "before",
      inclusive = TRUE
    )

  ovr <- filter(adrs, ANL01FL == "Y", ANL02FL == "Y")

  no_data_n <- event(dataset_name = "adsl", condition = TRUE,
                     set_values_to = exprs(AVALC = "N"), keep_source_vars = adsl_vars)
  no_data_missing <- event(dataset_name = "adsl", condition = TRUE,
                           set_values_to = exprs(AVALC = "MISSING"),
                           keep_source_vars = adsl_vars)
  src <- list(ovr = ovr, adsl = adsl)

  adrs <- adrs %>%
    # first progression (Y/N for every subject)
    derive_extreme_records(
      dataset_ref = adsl, dataset_add = adrs, by_vars = keys,
      filter_add = PARAMCD == "OVR" & AVALC == "PD" & ANL01FL == "Y",
      order = exprs(ADT, RSSEQ), mode = "first",
      exist_flag = AVALC, false_value = "N",
      set_values_to = exprs(PARAMCD = "PD", PARAM = "Disease Progression by Investigator",
                            !!!recist_labels, AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) %>%
    # response, no confirmation needed
    derive_extreme_event(
      by_vars = keys, order = exprs(event_nr, ADT), tmp_event_nr_var = event_nr,
      mode = "first", events = list(rsp_y, no_data_n), source_datasets = src,
      set_values_to = exprs(PARAMCD = "RSP",
                            PARAM = "Response by Investigator (confirmation not required)",
                            !!!recist_labels, AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) %>%
    # best overall response
    derive_extreme_event(
      by_vars = keys, order = exprs(event_nr, ADT), tmp_event_nr_var = event_nr,
      mode = "first", source_datasets = src,
      events = list(bor_cr, bor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, no_data_missing),
      set_values_to = exprs(PARAMCD = "BOR",
                            PARAM = "Best Overall Response by Investigator (confirmation not required)",
                            !!!recist_labels, AVAL = aval_resp(AVALC), ANL01FL = "Y")
    ) %>%
    # confirmed response and confirmed BOR
    derive_extreme_event(
      by_vars = keys, order = exprs(desc(AVALC), ADT, event_nr),
      tmp_event_nr_var = event_nr, mode = "first", source_datasets = src,
      events = list(crsp_y_cr, crsp_y_pr, no_data_n),
      set_values_to = exprs(PARAMCD = "CRSP", PARAM = "Confirmed Response by Investigator",
                            !!!recist_labels, AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) %>%
    derive_extreme_event(
      by_vars = keys, order = exprs(event_nr, ADT), tmp_event_nr_var = event_nr,
      mode = "first", source_datasets = src,
      events = list(cbor_cr, cbor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, no_data_missing),
      set_values_to = exprs(PARAMCD = "CBOR",
                            PARAM = "Best Confirmed Overall Response by Investigator",
                            !!!recist_labels, AVAL = aval_resp(AVALC), ANL01FL = "Y")
    )

  # death and last assessment are needed for the PFS/OS censoring rules
  adsldth <- select(adsl, !!!keys, DTHDT, !!!adsl_vars)
  adrs %>%
    derive_extreme_records(
      dataset_ref = adsldth, dataset_add = adsldth, by_vars = keys,
      filter_add = !is.na(DTHDT), exist_flag = AVALC, false_value = "N",
      set_values_to = exprs(PARAMCD = "DEATH", PARAM = "Death", PARCAT1 = "Reference Event",
                            ANL01FL = "Y", AVAL = yn_to_numeric(AVALC), ADT = DTHDT)
    ) %>%
    select(-DTHDT) %>%
    derive_extreme_records(
      dataset_ref = adsl, dataset_add = adrs, by_vars = keys,
      filter_add = PARAMCD == "OVR" & ANL01FL == "Y",
      order = exprs(ADT, RSSEQ), mode = "last",
      set_values_to = exprs(PARAMCD = "LSTA", PARAM = "Last Disease Assessment by Investigator",
                            !!!recist_labels, ANL01FL = "Y")
    ) %>%
    # ASEQ is needed later: ADTTE keeps it as the source record (SRCSEQ)
    derive_var_obs_number(by_vars = keys, order = exprs(PARAMCD, ADT, VISITNUM, RSSEQ),
                          check_type = "error")
}
