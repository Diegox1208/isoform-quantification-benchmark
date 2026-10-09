#!/usr/bin/env python3
#### created by Diego G. Salas
#
# Stage 02 prepares the reference for one chromosome or one region of it.
# It keeps every transcript of the selected gene types that lies whole inside the region and is long
# enough to simulate, and writes the genome, the matching GTF and the spliced transcript sequences.
# Each transcript gets the share of its 31-mers that no other isoform of the same gene contains.
# Coordinates are shifted so the region starts at 1.

import argparse
import gzip
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

##### Getopts #####
#-----------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--raw", required=True, help="folder with genome.fa.gz and genes.gtf.gz from stage 01")
parser.add_argument("--out", required=True, help="reference folder")
parser.add_argument("--chrom", default="chr22")
parser.add_argument("--start", type=int, default=None, help="1-based first base of the region, empty for the whole chromosome")
parser.add_argument("--end", type=int, default=None, help="1-based last base of the region")
parser.add_argument("--fasta-offset", type=int, default=0, help="bases missing before the first base of the raw FASTA")
parser.add_argument("--gene-types", default="protein_coding,lncRNA")
parser.add_argument("--min-length", type=int, default=400)
parser.add_argument("--kmer", type=int, default=31)
args = parser.parse_args()

raw = Path(args.raw)
out = Path(args.out)
for f in ("genome.fa.gz", "genes.gtf.gz"):
    if not (raw / f).is_file():
        print(f"/!\\ Error : {raw / f} does not exist, run stage 01 first", file=sys.stderr)
        sys.exit(1)
out.mkdir(parents=True, exist_ok=True)

##### Functions #####
#-----------------------------------------------------------
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path, chrom):
    opener = gzip.open if str(path).endswith(".gz") else open
    seq, keep = [], False
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0].split(":")[0] == chrom
            elif keep:
                seq.append(line.strip())
    return "".join(seq).upper()


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def attributes(field):
    pairs = (a.strip().split(" ", 1) for a in field.strip().rstrip(";").split(";") if a.strip())
    out = defaultdict(list)
    for k, v in pairs:
        out[k].append(v.strip('"'))
    return out


def kmers(seq, k):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


class Timer:
    def __init__(self, block):
        self.block = block

    def __enter__(self):
        self.t = time.time()

    def __exit__(self, *exc):
        print(f"{self.block}\tdone in {round(time.time() - self.t)}s", flush=True)

##### Data files #####
#-----------------------------------------------------------
with Timer("Genome"):
    chrom_seq = read_fasta(raw / "genome.fa.gz", args.chrom)
    if not chrom_seq:
        print(f"/!\\ Error : {args.chrom} not found in {raw / 'genome.fa.gz'}", file=sys.stderr)
        sys.exit(1)
    start = args.start or 1 + args.fasta_offset
    end = args.end or args.fasta_offset + len(chrom_seq)
    shift = start - 1
    genome = chrom_seq[start - 1 - args.fasta_offset:end - args.fasta_offset]
    if len(genome) != end - start + 1:
        print(f"/!\\ Error : region {start}-{end} falls outside the FASTA", file=sys.stderr)
        sys.exit(1)
    write_fasta(out / "genome.fa", [(args.chrom, genome)])
    print(f"{args.chrom}:{start}-{end}, {len(genome):,} bases, {genome.count('N'):,} N")

##### Transcripts #####
#-----------------------------------------------------------
with Timer("Transcripts"):
    gene_types = set(args.gene_types.split(","))
    transcripts = {}
    exons = defaultdict(list)
    lines = defaultdict(list)
    with gzip.open(raw / "genes.gtf.gz", "rt") as f:
        for line in f:
            if not line.startswith(args.chrom + "\t"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] not in ("transcript", "exon"):
                continue
            a = attributes(c[8])
            if a["gene_type"][0] not in gene_types:
                continue
            tid = a["transcript_id"][0]
            lines[tid].append(c)
            if c[2] == "transcript":
                transcripts[tid] = dict(gene_id=a["gene_id"][0], gene_name=a["gene_name"][0],
                                        transcript_type=a["transcript_type"][0],
                                        start=int(c[3]), end=int(c[4]), strand=c[6])
            else:
                exons[tid].append((int(c[3]), int(c[4])))

    records = []
    for tid, v in transcripts.items():
        if v["start"] < start or v["end"] > end:
            continue
        seq = "".join(genome[s - 1 - shift:e - shift] for s, e in sorted(exons[tid]))
        if v["strand"] == "-":
            seq = revcomp(seq)
        if len(seq) < args.min_length or "N" in seq:
            continue
        records.append(dict(transcript_id=tid, **v, length=len(seq), exons=len(exons[tid]), seq=seq))
    catalog = pd.DataFrame(records)
    catalog["isoforms"] = catalog.groupby("gene_id")["transcript_id"].transform("size")

    with open(out / "genes.gtf", "w") as f:
        for tid in catalog["transcript_id"]:
            for c in lines[tid]:
                c = c.copy()
                c[3], c[4] = str(int(c[3]) - shift), str(int(c[4]) - shift)
                f.write("\t".join(c) + "\n")
    print(f"{len(catalog):,} transcripts of {catalog['gene_id'].nunique():,} genes kept, "
          f"{(catalog['isoforms'] > 1).sum():,} of them share their gene with another isoform")

##### Unique sequence per isoform #####
#-----------------------------------------------------------
with Timer("Unique k-mers"):
    unique_share = {}
    for gene_id, g in catalog.groupby("gene_id"):
        sets = {t: kmers(s, args.kmer) for t, s in zip(g["transcript_id"], g["seq"])}
        counts = Counter(k for s in sets.values() for k in s)
        for t, s in sets.items():
            unique_share[t] = sum(counts[k] == 1 for k in s) / len(s) if s else 0.0
    catalog["unique_kmer_share"] = catalog["transcript_id"].map(unique_share).round(4)

##### Saving data #####
#-----------------------------------------------------------
write_fasta(out / "transcripts.fa", zip(catalog["transcript_id"], catalog["seq"]))
catalog.drop(columns="seq").to_csv(out / "catalog.tsv", sep="\t", index=False)
catalog[["gene_id", "transcript_id"]].to_csv(out / "tx2gene.tsv", sep="\t", index=False, header=False)
print(f"Saved genome.fa, genes.gtf, transcripts.fa, catalog.tsv and tx2gene.tsv in {out}")
