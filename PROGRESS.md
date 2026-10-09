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

### Test dataset

`bash tests/make_test_data.sh` → chr22 from 20 000 001 to 21 000 000, cut by the UCSC API, with the GENCODE lines that touch it, 0.3 MB in `tests/data/raw`. Two builds give the same bytes. The megabase keeps 80 transcripts of 37 genes, with up to 9 isoforms per gene, so it checks that every stage runs and nothing more.

### Result

Working run on the whole of chr22, nine samples of 1 000 000 fragments, 35 min 44 s on 6 threads.

<table>
<thead><tr><th>Item</th><th>Predicted</th><th>Outcome</th></tr></thead>
<tbody><tr><td>Salmon and kallisto per isoform at 2 × 100 bp</td><td>median absolute error under 10%</td><td>1.9% and 2.5%, confirmed</td></tr></tbody>
<tbody><tr><td>Isoform error from genes with 2 isoforms to genes with 5 or more</td><td>at least doubles</td><td>Salmon 1.9% to 11.7%, kallisto 2.3% to 12.3%, RSEM 4.2% to 12.4%, confirmed</td></tr></tbody>
<tbody><tr><td>All four methods per gene at 2 × 100 bp</td><td>median absolute error under 5%</td><td>Salmon 0.0%, kallisto 0.0%, RSEM 2.1%, featureCounts 0.2%, confirmed</td></tr></tbody>
<tbody><tr><td>Not predicted</td><td></td><td>featureCounts misses by more than 10% on 25.0% of genes and counts 92.7% of the fragments, while its median gene stays within 0.2%</td></tr></tbody>
<tbody><tr><td>Not predicted</td><td></td><td>RSEM counts 97.8% of the fragments at every read length, because STAR leaves soft-clipped alignments out of the transcriptome BAM</td></tr></tbody>
<tbody><tr><td>Not predicted</td><td></td><td>the sequence an isoform owns separates errors more sharply than the number of isoforms, Salmon 0.0% against 9.7%</td></tr></tbody>
<tbody><tr><td>Seed to seed spread</td><td></td><td>median isoform error moves 0.2 points or less, Spearman 0.004 or less</td></tr></tbody>
<tbody><tr><td><code>bash run.sh</code> on the test megabase</td><td>runs end to end</td><td>39 s on 2 threads</td></tr></tbody>
</table>

The three predictions held, which says they asked too little. The results worth keeping are the ones the prediction did not name, a gene-level median that hides a quarter of badly counted genes and a 2% loss that happens before RSEM sees the reads.

Errors met during the session, quoted from the logs with long paths cut, each fixed in its own commit.

- `No valid read to estimate the length distribution!` followed by `"rsem-run-em ..." failed!` on the test megabase, with polyester FASTA reads aligned by STAR → reads written as FASTQ with a constant quality.

### Costliest error

Escribí las lecturas en FASTA porque es lo que produce polyester, y STAR, Salmon y kallisto lo aceptaron sin quejarse. RSEM no. Sin calidades no encontró ninguna lectura válida para estimar la distribución de longitudes de fragmento y se detuvo con un mensaje que no menciona las calidades. Me costó unos 10 minutos entre revisar el BAM, leer los conteos de RSEM y probar a mano una muestra en FASTQ. La lección es que el formato de entrada se prueba contra la herramienta más exigente de la cadena y no contra la primera que lo acepta.

### Next step

Rerun RSEM after STAR in end-to-end mode to confirm that the 2.2% shortfall disappears. Then add a design with a flatter isoform split and with isoforms missing from the counting annotation, the error real data adds and this simulation leaves out.
