# The showcase trial as raw case report form (CRF) data, the way it would
# come out of an EDC system before any SDTM mapping: one data set per form,
# patient numbers instead of USUBJID, the EDC's own field names, collected
# text ("Female", "Germany") instead of CDISC terms, and dates as DD-Mon-YYYY text.
# sdtm.oak builds DM from these in trial/scripts/08_submission.R.

library(dplyr)

write_raw <- function(d, start = as.Date("2024-01-01")) {
  t <- pfs_os_truth(d)
  patnum <- sprintf("%04d", d$id)
  rand_dt <- start + d$rand_day
  # month names written out by hand so the dates don't depend on the locale
  crf_date <- function(x) sprintf("%02d-%s-%s", as.integer(format(x, "%d")),
                                  month.abb[as.integer(format(x, "%m"))], format(x, "%Y"))
  country_name <- c(USA = "United States", CAN = "Canada", DEU = "Germany", FRA = "France",
                    ESP = "Spain", GBR = "United Kingdom", CHN = "China", JPN = "Japan",
                    KOR = "South Korea")
  died <- t$os_event == 1

  list(
    # demographics form
    dm = tibble(PATNUM = patnum, AGE_YRS = d$age,
                GENDER = ifelse(d$sex == "F", "Female", "Male"),
                CNTRY = unname(country_name[d$country])),
    # randomisation form
    rand = tibble(PATNUM = patnum, RAND_DAT = crf_date(rand_dt),
                  TRT_ASSIGNED = ifelse(d$arm == 1, "Drug X", "Platinum chemotherapy")),
    # death form, only for patients who died
    death = tibble(PATNUM = patnum[died], DEATH_DAT = crf_date(rand_dt[died] + d$t_death[died]))
  )
}
