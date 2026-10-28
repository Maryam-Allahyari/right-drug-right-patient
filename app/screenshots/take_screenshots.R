# Screenshots of the teal app for the README. Start the app first
# (Rscript -e 'shiny::runApp("app", port = 8765)'), then from the repo root:
#   Rscript app/screenshots/take_screenshots.R
# chromote drives a headless Chrome or Edge; set CHROMOTE_CHROME to its path
# if it isn't found on its own.

library(chromote)

url <- "http://127.0.0.1:8765/"
tabs <- c(demographics = "demographics", km = "kaplan_meier_plot")

b <- ChromoteSession$new(width = 1500, height = 1000)
b$Page$navigate(url)
Sys.sleep(15)   # teal loads the data and renders the first module

for (nm in names(tabs)) {
  # open the module's tab, the same as clicking it in the menu
  b$Runtime$evaluate(sprintf(
    "document.querySelector('a[href=\"#teal-teal_modules-nav-%s-wrapper\"]').click()", tabs[[nm]]))
  Sys.sleep(25)   # long enough for the module to finish computing
  b$screenshot(file.path("app", "screenshots", paste0(nm, ".png")), show = FALSE)
}
b$close()
