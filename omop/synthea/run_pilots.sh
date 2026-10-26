#!/usr/bin/env bash
# Generates two small populations, one with lung cancer and one with breast
# cancer, so we can see which one has enough treatment and death records
# to build the cohort on.
#
# Before running:
#   - Java 17 and Synthea 3.3.0 (paths below, change them if yours differ).
#     Stay on 3.3.0: ETL-Synthea doesn't support newer versions yet.
#
# From anywhere:
#   bash omop/synthea/run_pilots.sh
#
# Output: Synthea CSVs in M:/Data/synthea/pilot_lung and pilot_breast
# (not tracked in git), plus a log for each. Takes about 20 minutes.
#
# A few things that surprised me:
#   - -p counts *living* patients only; dead ones come on top. Lung cancer
#     patients mostly die, so -p 25 still gave 204 patients.
#   - -a 45-90 keeps ages where these cancers actually happen, which saves
#     a lot of run time.
#   - the keep files (-k) throw away everyone without the cancer diagnosis.

set -e   # stop at the first error instead of carrying on

# where things are on my machine, change these to yours
JAVA="M:/Tools/jdk-17.0.20.1+1/bin/java"
JAR="M:/Tools/synthea/synthea-with-dependencies-3.3.0.jar"
KEEP="M:/Projects/right-drug-right-patient/omop/synthea"   # folder with the keep files
OUT="M:/Data/synthea"
mkdir -p "$OUT"

# run one pilot: run <cancer> <number of living patients>
run() {
  local name=$1 alive=$2
  echo "== $name: start $(date)"
  # -Xmx6g gives Java up to 6 GB of memory
  # -s 2026 fixes the random seed, so the same command gives the same patients
  # CSV only: that's what ETL-Synthea reads later, FHIR would just take up space
  "$JAVA" -Xmx6g -jar "$JAR" -s 2026 -p "$alive" -a 45-90 \
    -k "$KEEP/keep_$name.json" \
    --exporter.csv.export true --exporter.fhir.export false \
    --exporter.baseDirectory "$OUT/pilot_$name" > "$OUT/pilot_$name.log" 2>&1
  echo "== $name: done $(date)"
  tail -3 "$OUT/pilot_$name.log"   # last lines show total, alive and dead patients
}

run lung 25
run breast 120   # more living patients asked for, since far fewer breast cancer patients die
