# One trial taken through the whole pipeline, the way a real study would be
# reported: simulate -> SDTM -> ADaM (admiral/admiralonco) -> tables and
# figures with rtables/tern, numbered like a clinical study report.
# Scenario S4 (the effect grows with the biomarker), 800 patients.
#
# Before running:
#   install.packages(c("admiral", "admiralonco", "rtables", "tern",
#                      "ggsurvfit", "ggplot2", "survival", "dplyr", "tidyr"))
#
# From the repo root:
#   Rscript trial/scripts/05_showcase.R
#
# Output in trial/results/showcase/: tables as text, figures as png, the
# SDTM domains as csv, adtte_os.csv (OS per patient, used later to
# compare the control arm with the real-world cohort) and ADSL, ADRS and
# ADTTE as rds for the interactive app (app/app.R).

suppressPackageStartupMessages({
  library(dplyr); library(admiral); library(admiralonco)
  library(survival); library(splines); library(ggplot2); library(ggsurvfit)
  library(rtables); library(tern)
})
for (f in c("simulate", "write_sdtm", "derive_adsl", "derive_adrs", "derive_adtte"))
  source(file.path("trial/R", paste0(f, ".R")))

out <- "trial/results/showcase"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
arms <- c("Platinum chemotherapy", "Drug X")   # control first, it's the reference

d <- sim_trial(800, "S4", seed = 20261027)
sdtm <- write_sdtm(d)
for (nm in names(sdtm)) write.csv(sdtm[[nm]], file.path(out, paste0("sdtm_", nm, ".csv")), row.names = FALSE)

# factors fix the order of columns and rows in the tables
adsl <- derive_adsl(sdtm$dm, sdtm$ex, sdtm$ds, sdtm$sc, sdtm$rs, sdtm$ss) %>%
  mutate(ARM = factor(ARM, levels = arms), SEX = factor(SEX), ECOG = factor(ECOG),
         AGEGR1 = factor(AGEGR1, levels = c("<65", ">=65")),
         HISTOL = factor(HISTOL, levels = c("Non-squamous", "Squamous")),
         BMRKGR1 = factor(BMRKGR1, levels = c("Low", "High")),
         REGION1 = factor(REGION1, levels = c("North America", "Europe", "Asia")))
adrs  <- derive_adrs(sdtm$rs, adsl)
adtte_full <- derive_adtte(adsl, adrs)
adtte <- adtte_full %>%
  select(USUBJID, PARAMCD, AVAL, CNSR) %>%
  left_join(adsl, by = "USUBJID") %>%
  mutate(AVALM = AVAL / month, is_event = CNSR == 0)

os  <- filter(adtte, PARAMCD == "OS")
pfs <- filter(adtte, PARAMCD == "PFS")
write.csv(select(os, USUBJID, ARM, AGE, AGEGR1, SEX, ECOG, AVAL, AVALM, CNSR),
          file.path(out, "adtte_os.csv"), row.names = FALSE)

# the ADaM datasets as they are, for the teal app
saveRDS(adsl, file.path(out, "adsl.rds"))
saveRDS(adrs, file.path(out, "adrs.rds"))
saveRDS(adtte_full, file.path(out, "adtte.rds"))

# ---- Table 14.1.1 demographics and baseline characteristics ----
t_demo <- basic_table(show_colcounts = TRUE,
                      title = "Table 14.1.1 Demographics and Baseline Characteristics, ITT") %>%
  split_cols_by("ARM") %>% add_overall_col("All") %>%
  analyze_vars(c("AGE", "AGEGR1", "SEX", "REGION1", "ECOG", "HISTOL", "BMRK", "BMRKGR1")) %>%
  build_table(adsl)
writeLines(toString(t_demo), file.path(out, "t_14_1_1_demographics.txt"))

# time-to-event summary: events, median (KM), unstratified Cox HR
tte_table <- function(data, title, label) {
  basic_table(show_colcounts = TRUE, title = title) %>%
    split_cols_by("ARM", ref_group = arms[1]) %>%
    count_values("is_event", values = TRUE, .labels = c(count_fraction = "Patients with event")) %>%
    surv_time(vars = "AVALM", var_labels = label, is_event = "is_event",
              table_names = "time") %>%
    coxph_pairwise(vars = "AVALM", is_event = "is_event", var_labels = "Unstratified analysis",
                   table_names = "cox") %>%
    build_table(data)
}

# ---- Table 14.2.1 OS (primary) and Figure 14.2.1 KM plot ----
t_os <- tte_table(os, "Table 14.2.1 Overall Survival, ITT", "OS (months)")
writeLines(toString(t_os), file.path(out, "t_14_2_1_os.txt"))
print(t_os)

km <- survfit2(Surv(AVALM, is_event) ~ ARM, data = os) %>%
  ggsurvfit() + add_confidence_interval() + add_risktable() +
  labs(x = "Months since randomisation", y = "Overall survival",
       title = "Figure 14.2.1 Kaplan-Meier plot of OS, ITT")
ggsave(file.path(out, "f_14_2_1_km_os.png"), km, width = 8, height = 6, dpi = 150)

# ---- Table 14.2.2 PFS (secondary) ----
t_pfs <- tte_table(pfs, "Table 14.2.2 Progression-Free Survival, ITT", "PFS (months)")
writeLines(toString(t_pfs), file.path(out, "t_14_2_2_pfs.txt"))

# ---- Table 14.2.3 best overall response ----
bor <- filter(adrs, PARAMCD == "BOR") %>% select(USUBJID, AVALC) %>%
  left_join(select(adsl, USUBJID, ARM), by = "USUBJID") %>%
  mutate(AVALC = factor(AVALC, levels = c("CR", "PR", "SD", "PD", "NE", "MISSING")),
         RSP = AVALC %in% c("CR", "PR"))
t_bor <- basic_table(show_colcounts = TRUE,
                     title = "Table 14.2.3 Best Overall Response (investigator), ITT") %>%
  split_cols_by("ARM") %>%
  estimate_proportion(vars = "RSP", table_names = "orr") %>%
  analyze_vars("AVALC", .stats = "count_fraction", var_labels = c(AVALC = "Best overall response")) %>%
  build_table(bor)
writeLines(toString(t_bor), file.path(out, "t_14_2_3_response.txt"))

# ---- Figure 14.2.2 forest plot of OS by subgroup (tern) ----
# short arm labels here, the long one runs into "Drug X Better" in the header
df_sg <- extract_survival_subgroups(
  variables = list(tte = "AVALM", is_event = "is_event", arm = "ARM",
                   subgroups = c("SEX", "AGEGR1", "REGION1", "ECOG", "HISTOL", "BMRKGR1")),
  data = mutate(os, ARM = factor(ARM, levels = arms, labels = c("Chemo", "Drug X"))))
tbl_fp <- basic_table() %>%
  tabulate_survival_subgroups(df = df_sg, vars = c("n_tot", "n_tot_events", "hr", "ci"))
png(file.path(out, "f_14_2_2_forest_os.png"), width = 1800, height = 1300, res = 200)
print(g_forest(tbl_fp))
invisible(dev.off())

# ---- Figure 14.2.3 OS HR across the continuous biomarker, against the truth ----
# spline of log biomarker interacting with treatment, adjusted for the
# prognostic factors; the HR curve comes from the treatment terms only
os$lb <- log(os$BMRK)
kn <- quantile(os$lb, c(1/3, 2/3)); bd <- range(os$lb)
fit0 <- coxph(Surv(AVALM, is_event) ~ ARM + ns(lb, knots = kn, Boundary.knots = bd) +
                ECOG + AGEGR1 + HISTOL, data = os)
fit  <- coxph(Surv(AVALM, is_event) ~ ARM * ns(lb, knots = kn, Boundary.knots = bd) +
                ECOG + AGEGR1 + HISTOL, data = os)
grid <- exp(seq(log(quantile(os$BMRK, .05)), log(quantile(os$BMRK, .95)), length.out = 100))
b <- coef(fit); V <- vcov(fit)
idx <- c(grep("^ARMDrug X$", names(b)), grep("^ARMDrug X:", names(b)))
X <- cbind(1, ns(log(grid), knots = kn, Boundary.knots = bd))
est <- drop(X %*% b[idx]); se <- sqrt(rowSums((X %*% V[idx, idx]) * X))
curve <- data.frame(bmrk = grid, hr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
                    truth = exp(log(0.70) - 0.5 * log(grid / 5)))
p_bm <- ggplot(curve, aes(bmrk, hr)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_line(aes(y = truth), colour = "firebrick", linetype = 2) +
  geom_vline(xintercept = 5, linetype = 3) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Biomarker (log scale); dotted line = assay cut-off 5",
       y = "HR for OS, Drug X vs chemotherapy (log scale)",
       title = "Figure 14.2.3 Treatment effect on OS across biomarker levels",
       subtitle = "Blue: estimate adjusted for prognostic factors (95% CI). Red dashed: true HR") +
  theme_bw()
ggsave(file.path(out, "f_14_2_3_hr_by_biomarker.png"), p_bm, width = 8, height = 5, dpi = 150)
cat("biomarker x treatment interaction (spline) p =",
    signif(anova(fit0, fit)[2, "Pr(>|Chi|)"], 3), "\n")
