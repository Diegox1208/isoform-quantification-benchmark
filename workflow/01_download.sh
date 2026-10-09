#!/usr/bin/env bash
#### created by Diego G. Salas
#
# Stage 01 downloads the chromosome and the GENCODE annotation into ${refDir}/raw.
# A file already on disk is not downloaded again.
set -euo pipefail

#---- ##-- Configuration ----
source "${CONFIG:-tests/config.test.sh}"
DRY_RUN="${DRY_RUN:-0}"

##### Set variables
stage="01_download"
mkdir -p "${rawDir}"

#---- ##-- Loop over the two sources ----
for pair in "genome.fa.gz ${genomeUrl}" "genes.gtf.gz ${gtfUrl}"; do
  file="${pair%% *}"
  url="${pair#* }"
  out="${rawDir}/${file}"
  if [[ -s "${out}" ]]; then
    echo "[${stage}] ${out} exists, skipping" >&2
    continue
  fi

  ##### Build, log, run
  cmd="curl -fsSL --retry 3 -o ${out} '${url}'"
  echo "[${stage}] ${cmd}"
  echo -e "reference\t${stage}\t${cmd}\t${date_}" >> "${logFile}"
  if [[ "${DRY_RUN}" != "1" ]]; then
    eval "time ${cmd}"
  fi
done
