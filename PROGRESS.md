# PROGRESS

One entry per session, newest first. Structured fields in English, the reflection in Spanish.

---

## Session 1 · Simulated truth against four ways of quantifying isoforms

### Objective

Simulate paired-end fragments from every protein-coding and lncRNA isoform of chr22 with known counts, quantify them with Salmon, kallisto, RSEM on STAR and featureCounts at the gene level, and measure how far each method lands from the truth at 2 × 50, 2 × 100 and 2 × 150 bp, per isoform and per gene.

### Prediction

Written before any run.

1. At 2 × 100 bp, Salmon and kallisto stay under 10% median absolute error per isoform.
2. That error at least doubles from genes with 2 isoforms to genes with 5 or more.
3. At the gene level, all four methods stay under 5% median absolute error.
