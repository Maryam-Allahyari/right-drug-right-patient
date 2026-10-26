#!/usr/bin/env bash
# Pilot runs: small lung and breast cancer populations, to choose the disease.
# Synthea 3.3.0 (ETL-Synthea supports <= 3.3.0). Output goes to M:/Data/synthea (not in git).
#
# Notes
# - -p counts LIVING patients only; dead patients come on top. Lung cancer patients
#   mostly die, so a small -p still gives a few hundred patients.
# - -a 45-90 limits ages to where these cancers occur, which saves a lot of run time.
# - Keep files (-k) keep only patients with the cancer diagnosis.

set -e
JAVA="M:/Tools/jdk-17.0.20.1+1/bin/java"
JAR="M:/Tools/synthea/synthea-with-dependencies-3.3.0.jar"
KEEP="M:/Projects/right-drug-right-patient/omop/synthea"
OUT="M:/Data/synthea"
mkdir -p "$OUT"

run() {   # name, living patients
  local name=$1 alive=$2
  echo "== $name: start $(date)"
  "$JAVA" -Xmx6g -jar "$JAR" -s 2026 -p "$alive" -a 45-90 \
    -k "$KEEP/keep_$name.json" \
    --exporter.csv.export true --exporter.fhir.export false \
    --exporter.baseDirectory "$OUT/pilot_$name" > "$OUT/pilot_$name.log" 2>&1
  echo "== $name: done $(date)"
  tail -3 "$OUT/pilot_$name.log"
}

run lung 25
run breast 120
