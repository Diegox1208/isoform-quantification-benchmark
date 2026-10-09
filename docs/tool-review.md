# Tool review of Salmon, kallisto and RSEM for isoform quantification, benchmarked on reads simulated with polyester

Diego G. Salas · versions tested are Salmon 2.8.0, kallisto 0.52.0, RSEM 1.3.1 on STAR 2.7.11b, featureCounts 2.1.1 and polyester 1.38.0 on R 4.3.3 · dataset is GRCh38 chr22 with GENCODE release 46

## What it does

polyester draws fragments from a FASTA of transcripts with the exact number asked for each one, and writes the source transcript into every read name, so the truth can be read back from the output. Salmon and kallisto map each read pair straight to the transcript sequences and keep the set of transcripts it is compatible with. RSEM reads the alignments STAR writes against the transcriptome. All three then share every ambiguous fragment between its compatible isoforms with an expectation-maximisation step, which is the part that decides how a fragment shared by five isoforms is split. featureCounts counts fragments on genes from the genome alignments and leaves out any fragment that overlaps two genes or aligns to two places.

## Parameters that change the result

<table>
<thead><tr><th>Parameter</th><th>Tested values</th><th>Effect observed</th><th>Evidence</th></tr></thead>
<tbody><tr><td>Sequence the isoform owns</td><td>under 10% against 50% or more of its 31-mers</td><td>Salmon median error 9.7% against 0.0% at 2 × 100 bp, the same gradient for kallisto and RSEM</td><td><code>results/benchmark_strata.tsv</code></td></tr></tbody>
<tbody><tr><td>Isoforms in the gene</td><td>1, 2, 3-4, 5 or more</td><td>Salmon 0.0%, 1.9%, 5.6% and 11.7% at 2 × 100 bp</td><td><code>results/benchmark_strata.tsv</code></td></tr></tbody>
<tbody><tr><td>Read length</td><td>2 × 50, 2 × 100, 2 × 150 bp</td><td>kallisto median isoform error 3.4%, 2.5% and 2.4%, Salmon 2.3%, 1.9% and 1.9%</td><td><code>results/benchmark_summary.tsv</code></td></tr></tbody>
<tbody><tr><td>Counting level</td><td>featureCounts per gene against the isoform methods summed per gene</td><td>genes off by more than 10% at 25.0% against 1.2% to 2.9% at 2 × 100 bp</td><td><code>results/benchmark_summary.tsv</code></td></tr></tbody>
<tbody><tr><td>Alignments passed to RSEM</td><td>STAR transcriptome BAM with its default filter</td><td>97.8% of fragments counted at every read length</td><td>genome and transcriptome BAM files of <code>pe100_s1</code></td></tr></tbody>
<tbody><tr><td>Sample seed</td><td>1, 2, 3 per read length</td><td>standard deviation of the median isoform error 0.2 points or less, of the Spearman correlation 0.004 or less</td><td><code>results/benchmark_by_sample.tsv</code></td></tr></tbody>
</table>

## Where it fails

### Isoforms with little sequence of their own

When almost every 31-mer of an isoform also belongs to a sibling, only a few fragments tell them apart, and expectation-maximisation splits the rest by the abundances it has guessed. About a quarter of the isoforms with 20 true fragments or more end up off by more than 10% with all three methods.

### featureCounts on overlapping genes

Fragments that overlap two genes are left out as ambiguous and fragments that align twice as multimapping. The median gene stays within 0.2%, yet a quarter of genes miss by more than 10% and the Spearman correlation falls to 0.953. Longer reads overlap more genes, so the tail grows with read length.

### RSEM after STAR with soft-clipping

STAR soft-clips mismatched read ends by default, and the transcriptome BAM it writes for RSEM leaves those alignments out because RSEM cannot use them. In one sample at 2 × 100 bp, 19 439 of 997 386 aligned fragments never reach RSEM.

### RSEM on reads without qualities

With FASTA reads, STAR writes no base qualities and RSEM stops with `No valid read to estimate the length distribution!`. Salmon, kallisto and STAR accept the same files without complaint. Writing the reads as FASTQ with a constant quality solves it.

### kallisto on short reads

At 2 × 50 bp kallisto leaves 34.2% of isoforms off by more than 10%, against 26.8% for Salmon. The gap closes at 2 × 150 bp, 25.3% against 24.2%.

## Alternative

Salmon in alignment mode, run on the same STAR transcriptome BAM as RSEM, would separate the effect of mapping from the effect of the expectation-maximisation model. StringTie would test the case where the annotation is incomplete. Neither was run here.

## Source

- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 2017.
- Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology 2016.
- Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 2011.
- Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 2014.
- Dobin A, Davis CA, Schlesinger F et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 2013.
- Frazee AC, Jaffe AE, Langmead B, Leek JT. Polyester: simulating RNA-seq datasets with differential transcript expression. Bioinformatics 2015.
