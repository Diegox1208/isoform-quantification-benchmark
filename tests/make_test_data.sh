#!/usr/bin/env bash
#### created by Diego G. Salas
#
# Builds tests/data/raw, one megabase of chr22 in the same two formats stage 01 downloads.
# The sequence comes cut from the UCSC API, the GENCODE GTF is filtered here.
#   bash tests/make_test_data.sh
set -euo pipefail

#---- ##-- Configuration ----
source tests/config.test.sh
mkdir -p "${rawDir}"
api="https://api.genome.ucsc.edu/getData"
range="genome=hg38;chrom=${chrom};start=$((regionStart - 1));end=${regionEnd}"
gencode="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"

#---- ##-- Sequence, the region only, as genome.fa.gz ----
curl -fsSL "${api}/sequence?${range}" \
  | python -c "import json, sys, textwrap; d = json.load(sys.stdin)['dna']; print('>${chrom}'); print('\n'.join(textwrap.wrap(d.upper(), 60)))" \
  | gzip -n > "${rawDir}/genome.fa.gz"

#---- ##-- GENCODE lines that touch the region, as genes.gtf.gz ----
curl -fsSL "${gencode}" | gunzip -c \
  | awk -F'\t' -v c="${chrom}" -v s="${regionStart}" -v e="${regionEnd}" '$1 == c && $5 >= s && $4 <= e' \
  | gzip -n > "${rawDir}/genes.gtf.gz"

ls -l "${rawDir}"
du -sh "${rawDir}"
