# ADSL: one row per patient.
# Treatment, dates, disposition, death and baseline characteristics.
# Follows the admiral ADSL template, cut down to what this study uses.

library(admiral)
library(dplyr)

format_eosstt <- function(x) {
  case_when(
    x == "COMPLETED"      ~ "COMPLETED",
    x == "SCREEN FAILURE" ~ NA_character_,
    !is.na(x)             ~ "DISCONTINUED",
    TRUE                  ~ "ONGOING"
  )
}

region_of <- function(country) {
  case_when(country %in% c("USA", "CAN") ~ "North America",
            country %in% c("DEU", "FRA", "ESP", "GBR") ~ "Europe",
            country %in% c("CHN", "JPN", "KOR") ~ "Asia")
}

derive_adsl <- function(dm, ex, ds, sc, rs, ss) {
  dm <- convert_blanks_to_na(dm)
  ex <- convert_blanks_to_na(ex) %>%
    derive_vars_dtm(dtc = EXSTDTC, new_vars_prefix = "EXST") %>%
    derive_vars_dtm(dtc = EXENDTC, new_vars_prefix = "EXEN", time_imputation = "last")
  ds <- convert_blanks_to_na(ds) %>%
    derive_vars_dt(dtc = DSSTDTC, new_vars_prefix = "DSST")

  adsl <- dm %>%
    mutate(TRT01P = ARM, TRT01A = ACTARM) %>%
    # first and last dose of any study drug (both chemo drugs count)
    derive_vars_merged(
      dataset_add = ex, by_vars = exprs(STUDYID, USUBJID),
      filter_add = EXDOSE > 0 & !is.na(EXSTDTM),
      new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
      order = exprs(EXSTDTM, EXSEQ), mode = "first"
    ) %>%
    derive_vars_merged(
      dataset_add = ex, by_vars = exprs(STUDYID, USUBJID),
      filter_add = EXDOSE > 0 & !is.na(EXENDTM),
      new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
      order = exprs(EXENDTM, EXSEQ), mode = "last"
    ) %>%
    derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM)) %>%
    derive_var_trtdurd() %>%
    # disposition: end of study date and status
    derive_vars_merged(
      dataset_add = ds, by_vars = exprs(STUDYID, USUBJID),
      filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE",
      new_vars = exprs(EOSDT = DSSTDT)
    ) %>%
    derive_vars_merged(
      dataset_add = ds, by_vars = exprs(STUDYID, USUBJID),
      filter_add = DSCAT == "DISPOSITION EVENT",
      new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
      missing_values = exprs(EOSSTT = "ONGOING")
    ) %>%
    derive_vars_merged(
      dataset_add = ds, by_vars = exprs(STUDYID, USUBJID),
      filter_add = DSDECOD == "RANDOMIZED",
      new_vars = exprs(RANDDT = DSSTDT)
    ) %>%
    # death date; a partial date would be imputed to the 1st
    derive_vars_dt(dtc = DTHDTC, new_vars_prefix = "DTH",
                   highest_imputation = "M", date_imputation = "first") %>%
    derive_var_merged_exist_flag(
      dataset_add = ex, by_vars = exprs(STUDYID, USUBJID),
      new_var = SAFFL, condition = EXDOSE > 0,
      false_value = "N", missing_value = "N"
    ) %>%
    mutate(DTHFL = if_else(is.na(DTHFL), "N", DTHFL),
           ITTFL = "Y",                                  # everyone randomised
           AGEGR1 = if_else(AGE >= 65, ">=65", "<65"),
           REGION1 = region_of(COUNTRY))

  # baseline characteristics from SC, one column each
  base <- sc %>%
    select(USUBJID, SCTESTCD, SCSTRESC) %>%
    tidyr::pivot_wider(names_from = SCTESTCD, values_from = SCSTRESC) %>%
    mutate(BMRK = as.numeric(BMRK), BMRKGR1 = if_else(BMRK > 5, "High", "Low"))

  # last date known alive: latest of randomisation, any scan, last dose and
  # the last survival follow-up contact. OS is censored here for patients
  # who haven't died; without SS it would stop at the last scan.
  last_rs <- rs %>% group_by(USUBJID) %>%
    summarise(LASTRSDT = max(as.Date(RSDTC)), .groups = "drop")
  last_ss <- ss %>% filter(SSSTRESC == "ALIVE") %>% group_by(USUBJID) %>%
    summarise(LASTSSDT = max(as.Date(SSDTC)), .groups = "drop")

  adsl %>%
    left_join(base, by = "USUBJID") %>%
    left_join(last_rs, by = "USUBJID") %>%
    left_join(last_ss, by = "USUBJID") %>%
    mutate(LSTALVDT = pmax(RANDDT, LASTRSDT, TRTEDT, LASTSSDT, na.rm = TRUE)) %>%
    select(-LASTRSDT, -LASTSSDT)
}
