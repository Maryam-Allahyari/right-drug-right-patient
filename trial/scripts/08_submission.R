# The showcase trial's data made ready for a regulatory submission.
#
# SDTM: the DM domain built a second way, with sdtm.oak, from raw case report
# form (CRF) data, the way SDTM is built in a real study, and compared with
# the DM that write_sdtm.R writes directly.
#
# ADaM: a written specification (trial/specs/, one row per variable with
# label, type, length, format, origin, codelist and derivation, like a
# Pinnacle 21 spec) is read into metacore. metatools checks ADSL and ADTTE
# against it (no variable missing or extra, every value in its codelist),
# and xportr applies it and writes SAS transport files (.xpt v5), the
# format FDA and PMDA expect. Reading the files back checks nothing was lost.
#
# Before running:
#   - trial/scripts/05_showcase.R (adsl.rds and adtte.rds)
#   - install.packages(c("sdtm.oak", "metacore", "metatools", "xportr", "haven"))
#
# From the repo root:
#   Rscript trial/scripts/08_submission.R
#
# Output in trial/results/submission/: adsl.xpt, adtte.xpt, dm_oak.csv and
# checks.csv (what was checked and the result).

suppressPackageStartupMessages({
  library(dplyr); library(metacore); library(metatools); library(xportr); library(sdtm.oak)
})
source("trial/R/simulate.R")
source("trial/R/write_sdtm.R")
source("trial/R/write_raw.R")

out <- "trial/results/submission"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(what, ok) {
  checks[[length(checks) + 1]] <<- tibble(check = what, passed = ok)
  if (!ok) warning("check failed: ", what)
}

# ======================== SDTM DM with sdtm.oak ========================
d <- sim_trial(800, "S4", seed = 20261027)   # the showcase trial
raw <- write_raw(d)
ct <- read.csv("trial/specs/sdtm_ct.csv", colClasses = "character")

# every raw record gets the oak id variables, which link SDTM rows back to it
dm_raw   <- generate_oak_id_vars(raw$dm,   pat_var = "PATNUM", raw_src = "dm")
rand_raw <- generate_oak_id_vars(raw$rand, pat_var = "PATNUM", raw_src = "rand")
dth_raw  <- generate_oak_id_vars(raw$death, pat_var = "PATNUM", raw_src = "death")

# Each form is mapped on its own: sdtm.oak joins on the oak id variables,
# which belong to one raw data set. The mapped pieces are then joined by
# patient number, one DM row per patient.
dm_part <-
  assign_no_ct(raw_dat = dm_raw, raw_var = "PATNUM", tgt_var = "SUBJID", id_vars = oak_id_vars()) %>%
  assign_no_ct(raw_dat = dm_raw, raw_var = "AGE_YRS", tgt_var = "AGE", id_vars = oak_id_vars()) %>%
  # collected text ("Female", "Germany") mapped to CDISC terms ("F", "DEU")
  assign_ct(raw_dat = dm_raw, raw_var = "GENDER", tgt_var = "SEX",
            ct_spec = ct, ct_clst = "C66731", id_vars = oak_id_vars()) %>%
  assign_ct(raw_dat = dm_raw, raw_var = "CNTRY", tgt_var = "COUNTRY",
            ct_spec = ct, ct_clst = "C66786", id_vars = oak_id_vars())
rand_part <-
  assign_no_ct(raw_dat = rand_raw, raw_var = "TRT_ASSIGNED", tgt_var = "ARM", id_vars = oak_id_vars()) %>%
  assign_ct(raw_dat = rand_raw, raw_var = "TRT_ASSIGNED", tgt_var = "ARMCD",
            ct_spec = ct, ct_clst = "ARMCD", id_vars = oak_id_vars())
death_part <-
  assign_datetime(raw_dat = dth_raw, raw_var = "DEATH_DAT", tgt_var = "DTHDTC",
                  raw_fmt = "dd-mmm-yyyy", id_vars = oak_id_vars())

dm_oak <- dm_part %>%
  left_join(select(rand_part, patient_number, ARM, ARMCD), by = "patient_number") %>%
  left_join(select(death_part, patient_number, DTHDTC), by = "patient_number") %>%
  # study-level values, the same for every patient
  mutate(STUDYID = "NSCLC01", USUBJID = paste0("NSCLC01-", SUBJID), AGEU = "YEARS",
         ACTARM = ARM, ACTARMCD = ARMCD,
         DTHFL = if_else(!is.na(DTHDTC), "Y", NA_character_),
         DTHDTC = as.character(DTHDTC)) %>%
  select(STUDYID, USUBJID, SUBJID, ARM, ACTARM, ARMCD, ACTARMCD, AGE, AGEU, SEX,
         COUNTRY, DTHFL, DTHDTC) %>%
  arrange(USUBJID)
write.csv(dm_oak, file.path(out, "dm_oak.csv"), row.names = FALSE, na = "")

dm_direct <- write_sdtm(d)$dm %>% arrange(USUBJID)
same <- all.equal(as.data.frame(dm_oak), as.data.frame(dm_direct), check.attributes = FALSE)
print(same)
check("DM from raw CRF data with sdtm.oak is identical to DM from write_sdtm.R", isTRUE(same))

# ======================== ADaM: spec, checks, xpt ========================
# xportr counts "date" as a text type by default; SAS dates are numbers
# (days since 1960) with a date format, so "date" is moved to the numeric types
options(xportr.character_metadata_types = setdiff(getOption("xportr.character_metadata_types"), "date"),
        xportr.numeric_metadata_types = c(getOption("xportr.numeric_metadata_types"), "date"))
# ---- the spec as a metacore object ----
spec_ds  <- read.csv("trial/specs/adam_datasets.csv")
spec_var <- read.csv("trial/specs/adam_variables.csv", na.strings = "")
spec_cl  <- read.csv("trial/specs/adam_codelists.csv", colClasses = "character")

mc <- metacore(
  ds_spec = as_tibble(spec_ds),
  ds_vars = spec_var %>%
    transmute(dataset, variable, order, key_seq = as.integer(key_seq),
              mandatory = TRUE, core = NA_character_, supp_flag = FALSE),
  # variables used in both datasets are defined once ("common")
  var_spec = spec_var %>%
    group_by(variable) %>%
    summarise(label = first(label), type = first(type), length = first(length),
              format = first(format), common = n() > 1, .groups = "drop"),
  value_spec = spec_var %>%
    transmute(dataset, variable, origin, type, code_id = codelist, sig_dig = NA_integer_,
              derivation_id = if_else(is.na(derivation), NA_character_, paste0(dataset, ".", variable)),
              where = NA_character_),
  derivations = spec_var %>%
    filter(!is.na(derivation)) %>%
    transmute(derivation_id = paste0(dataset, ".", variable), derivation),
  codelist = spec_cl %>%
    group_by(code_id, name) %>%
    summarise(type = "code_decode", codes = list(tibble(code = code, decode = decode)),
              .groups = "drop"),
  supp = tibble(dataset = character(), variable = character(), idvar = character(),
                qeval = character())
)

# ---- check and apply the spec to one dataset, then write it ----
to_xpt <- function(data, ds) {
  ds_mc <- select_dataset(mc, ds, verbose = "silent")
  data <- data %>%
    mutate(across(where(is.factor), as.character)) %>%
    drop_unspec_vars(ds_mc)   # helper columns from the derivation that aren't in the spec

  # every specified variable is there, and every value is in its codelist
  check(paste(ds, "has exactly the variables in the spec"),
        isTRUE(tryCatch({ check_variables(data, ds_mc); TRUE }, error = function(e) FALSE)))
  check(paste(ds, "values are all in their codelists"),
        isTRUE(tryCatch({ check_ct_data(data, ds_mc, na_acceptable = TRUE); TRUE },
                        error = function(e) { message(conditionMessage(e)); FALSE })))
  check(paste(ds, "keys are unique"),
        isTRUE(tryCatch({ check_unique_keys(data, ds_mc); TRUE }, error = function(e) FALSE)))

  # type, length, label, order and format from the spec, then the file
  f <- file.path(out, paste0(tolower(ds), ".xpt"))
  ready <- data %>% order_cols(ds_mc) %>% sort_by_key(ds_mc)
  ready %>%
    xportr_type(ds_mc, domain = ds) %>%
    xportr_length(ds_mc, domain = ds) %>%
    xportr_label(ds_mc, domain = ds) %>%
    xportr_format(ds_mc, domain = ds) %>%
    xportr_df_label(ds_mc, domain = ds) %>%
    xportr_write(f)

  # read back: same values, labels in the file
  back <- haven::read_xpt(f)
  check(paste(ds, "xpt file has the spec's variables in the spec's order"),
        identical(names(back), filter(spec_var, dataset == ds)$variable))
  check(paste(ds, "xpt file has a label on every variable"),
        all(vapply(back, function(x) !is.null(attr(x, "label")), logical(1))))
  check(paste(ds, "xpt file has the same number of records"), nrow(back) == nrow(data))
  # values compared without their attributes; SAS stores missing text as ""
  plain <- function(x) { attributes(x) <- NULL; if (is.character(x)) x[x == ""] <- NA; x }
  check(paste(ds, "xpt file has the same values as the dataset"),
        all(vapply(names(back), function(v) isTRUE(all.equal(plain(back[[v]]), plain(ready[[v]]))),
                   logical(1))))
}

to_xpt(readRDS("trial/results/showcase/adsl.rds"), "ADSL")
to_xpt(readRDS("trial/results/showcase/adtte.rds"), "ADTTE")

checks <- bind_rows(checks)
print(checks)
write.csv(checks, file.path(out, "checks.csv"), row.names = FALSE)
