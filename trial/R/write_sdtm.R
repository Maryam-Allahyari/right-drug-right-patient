# Writes one simulated trial as SDTM domains: DM, DS, EX, RS, SC and SS.
# Only the variables the ADaM derivation needs; dates are ISO 8601 text,
# the way they arrive from a real study database.

library(dplyr)

write_sdtm <- function(d, study = "NSCLC01", start = as.Date("2024-01-01")) {
  t <- pfs_os_truth(d)
  usubjid <- sprintf("%s-%04d", study, d$id)
  rand_dt <- start + d$rand_day
  iso <- function(x) ifelse(is.na(x), NA_character_, format(x, "%Y-%m-%d"))
  death_obs <- t$os_event == 1
  arm_lbl <- ifelse(d$arm == 1, "Drug X", "Platinum chemotherapy")
  arm_cd  <- ifelse(d$arm == 1, "DRUGX", "CHEMO")

  dm <- tibble(
    STUDYID = study, USUBJID = usubjid, SUBJID = sprintf("%04d", d$id),
    ARM = arm_lbl, ACTARM = arm_lbl, ARMCD = arm_cd, ACTARMCD = arm_cd,
    AGE = d$age, AGEU = "YEARS", SEX = d$sex, COUNTRY = d$country,
    DTHFL = ifelse(death_obs, "Y", NA_character_),
    DTHDTC = iso(ifelse(death_obs, rand_dt + d$t_death, NA) |> as.Date(origin = "1970-01-01"))
  )

  # baseline characteristics; in a real study these come from several
  # domains (QS for ECOG, MI/LB for the biomarker...), one domain is enough here
  sc <- bind_rows(
    tibble(USUBJID = usubjid, SCTESTCD = "ECOG",   SCSTRESC = d$ecog),
    tibble(USUBJID = usubjid, SCTESTCD = "HISTOL", SCSTRESC = d$hist),
    tibble(USUBJID = usubjid, SCTESTCD = "BMRK",   SCSTRESC = as.character(d$bmrk))
  ) %>% mutate(STUDYID = study, .before = 1)

  # Drug X is one record; chemotherapy is cisplatin + paclitaxel, two records.
  # Doses are typical label doses, nothing downstream depends on them.
  ex_one <- function(i, trt, dose, seq) {
    tibble(STUDYID = study, USUBJID = usubjid[i], EXSEQ = seq, EXTRT = trt,
           EXDOSE = dose, EXDOSU = "mg/m2",
           EXSTDTC = iso(rand_dt[i]), EXENDTC = iso(rand_dt[i] + t$trt_end[i]))
  }
  drugx <- which(d$arm == 1); chemo <- which(d$arm == 0)
  ex <- bind_rows(
    ex_one(drugx, "DRUG X", 100, 1),
    ex_one(chemo, "CISPLATIN", 75, 1),
    ex_one(chemo, "PACLITAXEL", 175, 2)
  ) %>% arrange(USUBJID, EXSEQ)

  # randomisation, plus death or loss to follow-up as the disposition event
  lost <- !death_obs & d$dropped
  ds <- bind_rows(
    tibble(USUBJID = usubjid, DSCAT = "PROTOCOL MILESTONE", DSDECOD = "RANDOMIZED",
           DSSTDTC = iso(rand_dt)),
    tibble(USUBJID = usubjid[death_obs], DSCAT = "DISPOSITION EVENT", DSDECOD = "DEATH",
           DSSTDTC = iso(rand_dt[death_obs] + d$t_death[death_obs])),
    tibble(USUBJID = usubjid[lost], DSCAT = "DISPOSITION EVENT", DSDECOD = "LOST TO FOLLOW-UP",
           DSSTDTC = iso(rand_dt[lost] + d$fu_end[lost]))
  ) %>% mutate(STUDYID = study, .before = 1)

  # survival follow-up: the last contact for everyone still alive, at the
  # cut-off or when they were lost to follow-up
  alive <- !death_obs
  ss <- tibble(STUDYID = study, USUBJID = usubjid[alive], SSSEQ = 1,
               SSTESTCD = "SURVSTAT", SSTEST = "Survival Status", SSSTRESC = "ALIVE",
               SSDTC = iso(rand_dt[alive] + d$fu_end[alive]))

  # one investigator overall response per scan: SD, PR for responders from
  # the 2nd scan on, and PD at the scan that finds progression
  scans <- attr(d, "scans")
  rs <- lapply(which(d$n_scans > 0), function(i) {
    k <- seq_len(d$n_scans[i])
    resp <- ifelse(d$responder[i] & k >= 2, "PR", "SD")
    if (!is.na(d$pd_day[i])) resp[length(k)] <- "PD"
    tibble(USUBJID = usubjid[i], VISITNUM = k, VISIT = paste("WEEK", 6 * k),
           RSDTC = iso(rand_dt[i] + scans[i, k]), RSSTRESC = resp)
  }) %>% bind_rows() %>%
    mutate(STUDYID = study, RSTESTCD = "OVRLRESP", RSEVAL = "INVESTIGATOR") %>%
    group_by(USUBJID) %>% mutate(RSSEQ = row_number()) %>% ungroup()

  list(dm = dm, sc = sc, ex = ex, ds = ds, rs = rs, ss = ss)
}
