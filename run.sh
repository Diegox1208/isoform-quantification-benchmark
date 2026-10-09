#!/usr/bin/env bash
#### created by Diego G. Salas
#
# Runs the whole benchmark end to end.
#   bash run.sh                               -> one megabase of chr22, three samples, about a minute
#   CONFIG=config/config.sh bash run.sh       -> the whole of chr22, nine samples
#   DRY_RUN=1 bash run.sh                     -> print every command, run nothing
set -euo pipefail

#---- ##-- Configuration ----
export CONFIG="${CONFIG:-tests/config.test.sh}"
export DRY_RUN="${DRY_RUN:-0}"
source "${CONFIG}"

##### Build, log, run one stage
run_stage() {
  local stage="$1"; shift
  local cmd="$*"
  echo "[${stage}] ${cmd}"
  echo -e "all\t${stage}\t${cmd}\t${date_}" >> "${logFile}"
  if [[ "${DRY_RUN}" != "1" ]]; then
    eval "time ${cmd}"
  fi
}

regionArgs=""
if [[ -n "${regionStart}" ]]; then
  regionArgs="--start ${regionStart} --end ${regionEnd} --fasta-offset ${fastaOffset}"
fi

#---- ##-- Stages ----
bash workflow/01_download.sh

run_stage 02_prepare python workflow/02_prepare.py --raw "${rawDir}" --out "${refDir}" --chrom "${chrom}" ${regionArgs} \
  --gene-types "${geneTypes}" --min-length "${minTranscriptLength}"

run_stage 03_design python workflow/03_design.py --catalog "${refDir}/catalog.tsv" --out "${simDir}/design.tsv" \
  --fragments "${fragmentsPerSample}" --expressed "${expressedGeneFraction}" --alpha "${isoformAlpha}" --seed "${truthSeed}"

run_stage 04_simulate Rscript workflow/04_simulate.R --fasta "${refDir}/transcripts.fa" --design "${simDir}/design.tsv" \
  --info "${infoFile}" --out "${simDir}" --frag-mean "${fragmentMean}" --frag-sd "${fragmentSd}" --error-rate "${errorRate}" \
  --log "${logFile}" --stamp "${date_}"

bash workflow/05_index_align.sh
bash workflow/06_quantify.sh

run_stage 07_benchmark python workflow/07_benchmark.py --info "${infoFile}" --ref "${refDir}" --sim "${simDir}" \
  --quant "${quantDir}" --prefix "${resultsDir}/benchmark" --min-truth "${minTruth}"

echo "done, results in ${resultsDir}, log in ${logFile}"
