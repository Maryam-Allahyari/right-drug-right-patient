#!/usr/bin/env bash
# Generates the full lung cancer population that the cohort is built from.
# Lung was picked from the pilots (see omop/README.md); the settings are the
# same as the pilot, just more patients.
#
# Before running:
#   - Java 17 and Synthea 3.3.0 (paths below, change them if yours differ)
#   - about 6 GB free disk space
#
# From anywhere:
#   bash omop/synthea/run_full.sh
#
# Output: Synthea CSVs in M:/Data/synthea/full_lung (not tracked in git),
# plus a log in M:/Data/synthea/full_lung.log. Took about 25 minutes.
#
# -p 120 asks for 120 *living* patients. Most lung cancer patients in Synthea
# die, so the total is much bigger: 835 patients (714 NSCLC, 121 SCLC).
# Some living-patient slots give up after 1,000 failed attempts; that's
# expected for a rare disease and just means fewer living patients.

set -e   # stop at the first error instead of carrying on

# where things are on my machine, change these to yours
JAVA="M:/Tools/jdk-17.0.20.1+1/bin/java"
JAR="M:/Tools/synthea/synthea-with-dependencies-3.3.0.jar"
KEEP="M:/Projects/right-drug-right-patient/omop/synthea"   # folder with the keep files
OUT="M:/Data/synthea"
mkdir -p "$OUT"

echo "start $(date)"
# same options as the pilot (see run_pilots.sh), only -p is bigger and
# Java gets 8 GB instead of 6 because the run is larger
"$JAVA" -Xmx8g -jar "$JAR" -s 2026 -p 120 -a 45-90 \
  -k "$KEEP/keep_lung.json" \
  --exporter.csv.export true --exporter.fhir.export false \
  --exporter.baseDirectory "$OUT/full_lung" > "$OUT/full_lung.log" 2>&1
echo "done $(date)"
tail -3 "$OUT/full_lung.log"   # last lines show total, alive and dead patients
du -sh "$OUT/full_lung"         # how much disk space the output took
