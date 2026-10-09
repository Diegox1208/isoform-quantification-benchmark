#!/usr/bin/env bash
#### created by Diego G. Salas
#
# Stage 05 builds every index once, Salmon and kallisto on the transcripts, STAR and RSEM on the
# genome with the same GTF, then aligns each sample with STAR to the genome and to the transcriptome.
set -euo pipefail

#---- ##-- Configuration ----
source "${CONFIG:-tests/config.test.sh}"
DRY_RUN="${DRY_RUN:-0}"

##### Set variables
stage="05_index_align"
idxDir="${refDir}/index"
mkdir -p "${idxDir}"

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

#---- ##-- Indexes ----
if [[ ! -s "${idxDir}/salmon/versionInfo.json" ]]; then
  run index "salmon index -t ${refDir}/transcripts.fa -i ${idxDir}/salmon -k 31 -p ${threads}"
fi
if [[ ! -s "${idxDir}/kallisto.idx" ]]; then
  run index "kallisto index -i ${idxDir}/kallisto.idx -t ${threads} ${refDir}/transcripts.fa"
fi
if [[ ! -s "${idxDir}/star/SA" ]]; then
  mkdir -p "${idxDir}/star"
  run index "STAR --runMode genomeGenerate --runThreadN ${threads} --genomeDir ${idxDir}/star \
    --genomeFastaFiles ${refDir}/genome.fa --sjdbGTFfile ${refDir}/genes.gtf --sjdbOverhang 149 \
    --genomeSAindexNbases ${genomeSAindexNbases} --outFileNamePrefix ${idxDir}/star/ \
    --outTmpDir ${tmpDir}/STARtmp_index_$$"
  rm -rf "${tmpDir}/STARtmp_index_$$"
fi
if [[ ! -s "${idxDir}/rsem/ref.grp" ]]; then
  mkdir -p "${idxDir}/rsem"
  run index "rsem-prepare-reference --gtf ${refDir}/genes.gtf ${refDir}/genome.fa ${idxDir}/rsem/ref"
fi

#---- ##-- Loop over samples ----
while IFS=$'\t' read -r sample name readLength seed source; do
  [[ "${sample}" == "sample" || -z "${sample}" ]] && continue

  r1="${simDir}/${name}/sample_01_1.fasta.gz"
  r2="${simDir}/${name}/sample_01_2.fasta.gz"
  if [[ ! -f "${r1}" && "${DRY_RUN}" != "1" ]]; then
    echo "[${stage}] missing ${r1}, skipping ${name}" >&2
    continue
  fi

  ##### Name every output before running anything
  outDir="${alignDir}/${name}"
  prefix="${outDir}/${name}_"
  if [[ -s "${prefix}Aligned.toTranscriptome.out.bam" ]]; then
    echo "[${stage}] ${prefix}Aligned.toTranscriptome.out.bam exists, skipping ${name}" >&2
    continue
  fi
  mkdir -p "${outDir}"

  run "${name}" "STAR --runThreadN ${threads} --genomeDir ${idxDir}/star \
    --readFilesIn ${r1} ${r2} --readFilesCommand zcat \
    --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
    --outFileNamePrefix ${prefix} --outTmpDir ${tmpDir}/STARtmp_${name}_$$"
  rm -rf "${tmpDir}/STARtmp_${name}_$$"
done < "${infoFile}"
