#!/usr/bin/env bash
#### created by Diego G. Salas
#
# Stage 06 quantifies every sample with four methods.
# Salmon and kallisto map the reads straight to the transcripts, RSEM reads the STAR transcriptome BAM,
# and featureCounts counts the STAR genome BAM at the gene level only, as the baseline.
set -euo pipefail

#---- ##-- Configuration ----
source "${CONFIG:-tests/config.test.sh}"
DRY_RUN="${DRY_RUN:-0}"

##### Set variables
stage="06_quantify"
idxDir="${refDir}/index"

##### Build, log, run one command
run() {
  local name="$1"; shift
  local cmd="$*"
  echo "[${stage}] ${cmd}"
  echo -e "${name}\t${stage}\t${cmd}\t${date_}" >> "${logFile}"
  if [[ "${DRY_RUN}" != "1" ]]; then
    eval "time ${cmd}"
  fi
}

#---- ##-- Salmon, kallisto and RSEM, one folder per sample and method ----
bams=()
while IFS=$'\t' read -r sample name readLength seed source; do
  [[ "${sample}" == "sample" || -z "${sample}" ]] && continue
  r1="${simDir}/${name}/sample_01_1.fasta.gz"
  r2="${simDir}/${name}/sample_01_2.fasta.gz"
  prefix="${alignDir}/${name}/${name}_"
  if [[ ! -s "${prefix}Aligned.toTranscriptome.out.bam" && "${DRY_RUN}" != "1" ]]; then
    echo "[${stage}] missing ${prefix}Aligned.toTranscriptome.out.bam, skipping ${name}" >&2
    continue
  fi
  bams+=("${prefix}Aligned.out.bam")
  out="${quantDir}/${name}"
  mkdir -p "${out}"

  if [[ ! -s "${out}/salmon/quant.sf" ]]; then
    run "${name}" "salmon quant -i ${idxDir}/salmon -l A -1 ${r1} -2 ${r2} -p ${threads} \
      -o ${out}/salmon --quiet"
  fi
  if [[ ! -s "${out}/kallisto/abundance.tsv" ]]; then
    run "${name}" "kallisto quant -i ${idxDir}/kallisto.idx -o ${out}/kallisto -t ${threads} ${r1} ${r2}"
  fi
  if [[ ! -s "${out}/rsem.isoforms.results" ]]; then
    run "${name}" "rsem-calculate-expression --alignments --paired-end --no-bam-output -p ${threads} \
      --temporary-folder ${tmpDir}/rsem_${name}_$$ \
      ${prefix}Aligned.toTranscriptome.out.bam ${idxDir}/rsem/ref ${out}/rsem"
  fi
done < "${infoFile}"

#---- ##-- featureCounts, gene level, all samples in one table ----
out="${quantDir}/featureCounts_genes.txt"
if [[ -s "${out}" ]]; then
  echo "[${stage}] ${out} exists, skipping featureCounts" >&2
else
  run featureCounts "featureCounts -T ${threads} -p --countReadPairs -s 0 -t exon -g gene_id \
    -a ${refDir}/genes.gtf -o ${out} ${bams[*]}"
fi
